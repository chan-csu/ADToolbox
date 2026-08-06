"""Gene-marker evidence to e-ADM COD-group allocation.

The classifier is intentionally independent of the annotation program.  Its input
is a tall table of genome/marker hits; aliases and biological rules live in a
versioned JSON catalog so they can be reviewed and extended without changing
Python code.
"""

from __future__ import annotations

import json
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import polars as pl


DEFAULT_MARKER_CATALOG = Path(__file__).with_name("pkg_data") / "gene_marker_catalog.json"


class MarkerCatalogError(ValueError):
    """Raised when a marker catalog is internally inconsistent."""


@dataclass(frozen=True)
class MarkerPanel:
    id: str
    marker_weights: Mapping[str, float]
    required_all: tuple[str, ...]
    required_any: tuple[tuple[str, ...], ...]
    minimum_markers: int
    minimum_score: float
    strict: bool


@dataclass(frozen=True)
class MarkerGroup:
    id: str
    description: str
    panels: tuple[MarkerPanel, ...]


@dataclass(frozen=True)
class MarkerCatalog:
    schema_version: int
    catalog_version: str
    cod_groups: tuple[str, ...]
    markers: Mapping[str, Mapping[str, Any]]
    groups: Mapping[str, MarkerGroup]
    alias_lookup: Mapping[str, str]

    @classmethod
    def from_json(cls, path: str | Path = DEFAULT_MARKER_CATALOG) -> "MarkerCatalog":
        with Path(path).open(encoding="utf-8") as handle:
            payload = json.load(handle)
        return cls.from_dict(payload)

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "MarkerCatalog":
        if payload.get("schema_version") != 1:
            raise MarkerCatalogError("Only marker catalog schema_version 1 is supported")

        cod_groups = tuple(payload.get("cod_groups", ()))
        if not cod_groups or len(set(cod_groups)) != len(cod_groups):
            raise MarkerCatalogError("cod_groups must be a non-empty list of unique names")

        marker_data = payload.get("markers", {})
        if not isinstance(marker_data, dict) or not marker_data:
            raise MarkerCatalogError("markers must be a non-empty object")

        alias_lookup: dict[str, str] = {}
        for marker_id, definition in marker_data.items():
            aliases = [marker_id, *definition.get("aliases", [])]
            for alias in aliases:
                key = str(alias).strip().casefold()
                if not key:
                    raise MarkerCatalogError(f"Marker {marker_id!r} has an empty alias")
                if key in alias_lookup and alias_lookup[key] != marker_id:
                    raise MarkerCatalogError(
                        f"Alias {alias!r} is assigned to both {alias_lookup[key]!r} and {marker_id!r}"
                    )
                alias_lookup[key] = marker_id

        raw_groups = payload.get("groups", {})
        unknown_groups = set(raw_groups) - set(cod_groups)
        if unknown_groups:
            raise MarkerCatalogError(f"Rules reference unknown COD groups: {sorted(unknown_groups)}")

        groups: dict[str, MarkerGroup] = {}
        for group_id, group_data in raw_groups.items():
            panels: list[MarkerPanel] = []
            panel_ids: set[str] = set()
            for raw_panel in group_data.get("panels", []):
                panel_id = str(raw_panel.get("id", "")).strip()
                if not panel_id or panel_id in panel_ids:
                    raise MarkerCatalogError(f"COD group {group_id} has an empty or duplicate panel id")
                panel_ids.add(panel_id)

                weights = {str(k): float(v) for k, v in raw_panel.get("marker_weights", {}).items()}
                if not weights or any(value <= 0 for value in weights.values()):
                    raise MarkerCatalogError(f"Panel {panel_id} needs positive marker_weights")
                unknown_markers = set(weights) - set(marker_data)
                if unknown_markers:
                    raise MarkerCatalogError(
                        f"Panel {panel_id} references unknown markers: {sorted(unknown_markers)}"
                    )

                required_all = tuple(raw_panel.get("required_all", ()))
                required_any = tuple(tuple(clause) for clause in raw_panel.get("required_any", ()))
                requirement_markers = set(required_all)
                for clause in required_any:
                    if not clause:
                        raise MarkerCatalogError(f"Panel {panel_id} has an empty required_any clause")
                    requirement_markers.update(clause)
                missing_weights = requirement_markers - set(weights)
                if missing_weights:
                    raise MarkerCatalogError(
                        f"Panel {panel_id} requirements lack weights: {sorted(missing_weights)}"
                    )

                minimum_markers = int(raw_panel.get("minimum_markers", 1))
                minimum_score = float(raw_panel.get("minimum_score", 0.0))
                strict = bool(raw_panel.get("strict", False))
                if not 1 <= minimum_markers <= len(weights):
                    raise MarkerCatalogError(f"Panel {panel_id} has invalid minimum_markers")
                if not 0 <= minimum_score <= 1:
                    raise MarkerCatalogError(f"Panel {panel_id} has invalid minimum_score")

                panels.append(
                    MarkerPanel(
                        id=panel_id,
                        marker_weights=weights,
                        required_all=required_all,
                        required_any=required_any,
                        minimum_markers=minimum_markers,
                        minimum_score=minimum_score,
                        strict=strict,
                    )
                )
            if not panels:
                raise MarkerCatalogError(f"Configured COD group {group_id} has no panels")
            groups[group_id] = MarkerGroup(
                id=group_id,
                description=str(group_data.get("description", "")),
                panels=tuple(panels),
            )

        return cls(
            schema_version=1,
            catalog_version=str(payload.get("catalog_version", "unversioned")),
            cod_groups=cod_groups,
            markers=marker_data,
            groups=groups,
            alias_lookup=alias_lookup,
        )


def read_table(path: str | Path) -> pl.DataFrame:
    """Read a CSV or tab-delimited marker/abundance table."""
    path = Path(path)
    separator = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pl.read_csv(path, separator=separator)


def alignment_marker_hits(
    alignment: pl.DataFrame | str | Path,
    catalog: MarkerCatalog,
    *,
    entity_id: str,
    min_bits: float | None = None,
    max_evalue: float | None = None,
) -> pl.DataFrame:
    """Extract canonical marker hits from an MMseqs format-mode 4 table.

    Target FASTA headers may be ``accession|marker_id`` or contain the marker
    alias in any pipe-delimited field. Multiple reference sequences for the
    same family are collapsed to one query/marker observation.
    """
    frame = read_table(alignment) if isinstance(alignment, (str, Path)) else alignment
    required = {"query", "target"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"alignment is missing required columns: {sorted(missing)}")
    if frame.height == 0:
        return pl.DataFrame(schema={"genome_id": pl.Utf8, "gene_id": pl.Utf8, "marker_id": pl.Utf8})
    if min_bits is not None:
        if "bits" not in frame.columns:
            raise ValueError("alignment is missing required column: bits")
        frame = frame.filter(pl.col("bits").cast(pl.Float64, strict=False) >= min_bits)
    if max_evalue is not None:
        if "evalue" not in frame.columns:
            raise ValueError("alignment is missing required column: evalue")
        frame = frame.filter(pl.col("evalue").cast(pl.Float64, strict=False) <= max_evalue)

    rows: list[dict[str, Any]] = []
    optional = [column for column in ("bits", "evalue", "fident", "qcov", "tcov") if column in frame.columns]
    for row in frame.select("query", "target", *optional).to_dicts():
        target_tokens = [token.strip().casefold() for token in str(row["target"]).split("|")]
        canonical = next((catalog.alias_lookup[token] for token in target_tokens if token in catalog.alias_lookup), None)
        if canonical is None:
            continue
        record = {
            "genome_id": str(entity_id),
            "gene_id": str(row["query"]),
            "marker_id": canonical,
        }
        record.update({column: row[column] for column in optional})
        rows.append(record)
    schema = {"genome_id": pl.Utf8, "gene_id": pl.Utf8, "marker_id": pl.Utf8}
    if not rows:
        return pl.DataFrame(schema=schema)
    return pl.DataFrame(rows).unique(["genome_id", "gene_id", "marker_id"], keep="first")


def marker_ids_from_fasta(path: str | Path, catalog: MarkerCatalog) -> set[str]:
    """Return canonical markers represented in pipe-delimited FASTA headers."""
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    found: set[str] = set()
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip().split(None, 1)[0]
            tokens = [token.strip().casefold() for token in header.split("|")]
            found.update(catalog.alias_lookup[token] for token in tokens if token in catalog.alias_lookup)
    return found


def _resolve_profile_marker(value: str, catalog: MarkerCatalog) -> str | None:
    key = str(value).strip().casefold()
    direct = catalog.alias_lookup.get(key)
    if direct is not None:
        return direct
    return catalog.alias_lookup.get(key.split("__", 1)[0])


def hmmer_marker_hits(
    domtblout: str | Path,
    catalog: MarkerCatalog,
    *,
    entity_id: str,
    max_evalue: float = 1e-15,
    min_coverage: float = 0.35,
    score_cutoffs: Mapping[str, float | Mapping[str, Any]] | None = None,
) -> pl.DataFrame:
    """Parse ``hmmsearch --domtblout`` into canonical marker hits."""
    rows: list[dict[str, Any]] = []
    score_cutoffs = score_cutoffs or {}
    with Path(domtblout).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                continue
            gene_id = fields[0]
            profile_name = fields[3]
            profile_accession = fields[4]
            marker_id = _resolve_profile_marker(profile_name, catalog) or _resolve_profile_marker(profile_accession, catalog)
            if marker_id is None:
                continue
            query_length = float(fields[5])
            full_score = float(fields[7])
            domain_score = float(fields[13])
            independent_evalue = float(fields[12])
            hmm_from = float(fields[15])
            hmm_to = float(fields[16])
            coverage = (hmm_to - hmm_from + 1.0) / query_length if query_length > 0 else 0.0
            cutoff = score_cutoffs.get(profile_name, score_cutoffs.get(marker_id, float("-inf")))
            if isinstance(cutoff, Mapping):
                marker_cutoff = float(cutoff.get("score_threshold", float("-inf")))
                score_type = str(cutoff.get("score_type", "full")).strip().casefold()
            else:
                marker_cutoff = float(cutoff)
                score_type = "full"
            if score_type not in {"full", "domain"}:
                raise ValueError(f"Unsupported HMM score_type {score_type!r} for profile {profile_name!r}")
            selected_score = domain_score if score_type == "domain" else full_score
            if independent_evalue > max_evalue or coverage < min_coverage or selected_score < marker_cutoff:
                continue
            rows.append(
                {
                    "genome_id": str(entity_id),
                    "gene_id": gene_id,
                    "marker_id": marker_id,
                    "profile_id": profile_name,
                    "bits": selected_score,
                    "evalue": independent_evalue,
                    "coverage": coverage,
                }
            )
    schema = {
        "genome_id": pl.Utf8,
        "gene_id": pl.Utf8,
        "marker_id": pl.Utf8,
        "profile_id": pl.Utf8,
        "bits": pl.Float64,
        "evalue": pl.Float64,
        "coverage": pl.Float64,
    }
    if not rows:
        return pl.DataFrame(schema=schema)
    return pl.DataFrame(rows, schema=schema).sort("bits", descending=True).unique(
        ["genome_id", "gene_id", "marker_id"], keep="first"
    )


def read_hmm_score_cutoffs(path: str | Path | None) -> dict[str, dict[str, Any]]:
    """Read adaptive HMM cutoffs, including KOfam full/domain score type."""
    if path is None:
        return {}
    table = read_table(path)
    missing = {"profile_id", "score_threshold"} - set(table.columns)
    if missing:
        raise ValueError(f"HMM cutoff table is missing required columns: {sorted(missing)}")
    columns = ["profile_id", "score_threshold"]
    has_score_type = "score_type" in table.columns
    if has_score_type:
        columns.append("score_type")
    cutoffs: dict[str, dict[str, Any]] = {}
    for row in table.select(columns).drop_nulls(subset=["profile_id", "score_threshold"]).to_dicts():
        score_type = str(row.get("score_type") or "full").strip().casefold()
        if score_type not in {"full", "domain"}:
            raise ValueError(f"Unsupported HMM score_type {score_type!r} for profile {row['profile_id']!r}")
        cutoffs[str(row["profile_id"])] = {
            "score_threshold": float(row["score_threshold"]),
            "score_type": score_type,
        }
    return cutoffs


def build_marker_protein_db(
    manifest_path: str | Path,
    output_path: str | Path,
    catalog: MarkerCatalog,
) -> dict[str, int]:
    """Combine curated per-family FASTAs into an MMseqs marker target FASTA.

    The CSV/TSV manifest requires ``marker_id`` and ``source_fasta``. Relative
    source paths are resolved beside the manifest. Output headers are rewritten
    as ``source_accession|canonical_marker_id`` while descriptions are retained.
    """
    manifest_path = Path(manifest_path)
    manifest = read_table(manifest_path)
    missing = {"marker_id", "source_fasta"} - set(manifest.columns)
    if missing:
        raise ValueError(f"marker DB manifest is missing required columns: {sorted(missing)}")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    with output_path.open("w", encoding="utf-8") as output:
        for row in manifest.select("marker_id", "source_fasta").to_dicts():
            canonical = catalog.alias_lookup.get(str(row["marker_id"]).strip().casefold())
            if canonical is None:
                raise ValueError(f"Unknown catalog marker in manifest: {row['marker_id']!r}")
            source = Path(str(row["source_fasta"])).expanduser()
            if not source.is_absolute():
                source = manifest_path.parent / source
            if not source.exists():
                raise FileNotFoundError(source)
            opener = gzip.open if source.suffix.lower() == ".gz" else open
            sequence_count = 0
            with opener(source, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(">"):
                        header = line[1:].strip()
                        accession, _, description = header.partition(" ")
                        output.write(f">{accession}|{canonical}")
                        if description:
                            output.write(f" {description}")
                        output.write("\n")
                        sequence_count += 1
                    else:
                        output.write(line)
            if sequence_count == 0:
                raise ValueError(f"Reference FASTA contains no records: {source}")
            counts[canonical] = counts.get(canonical, 0) + sequence_count
    return counts


def build_marker_hmm_db(
    manifest_path: str | Path,
    output_path: str | Path,
    catalog: MarkerCatalog,
    *,
    cutoff_output: str | Path | None = None,
) -> dict[str, int]:
    """Extract and relabel selected HMM profiles from source HMM collections.

    The manifest requires ``marker_id``, ``source_hmm``, and ``profile_id``.
    An optional ``score_threshold`` is preserved in a companion cutoff table.
    Source HMM paths may point to dbCAN, KOfam, Pfam, or custom HMMER files.
    """
    manifest_path = Path(manifest_path)
    manifest = read_table(manifest_path)
    required = {"marker_id", "source_hmm", "profile_id"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"marker HMM manifest is missing required columns: {sorted(missing)}")

    requests_by_source: dict[Path, dict[str, tuple[str, float | None, str]]] = {}
    for row in manifest.to_dicts():
        canonical = catalog.alias_lookup.get(str(row["marker_id"]).strip().casefold())
        if canonical is None:
            raise ValueError(f"Unknown catalog marker in HMM manifest: {row['marker_id']!r}")
        source = Path(str(row["source_hmm"])).expanduser()
        if not source.is_absolute():
            source = manifest_path.parent / source
        if not source.exists():
            raise FileNotFoundError(source)
        profile_id = str(row["profile_id"]).strip()
        if not profile_id:
            raise ValueError("Every HMM manifest row needs a profile_id")
        threshold = row.get("score_threshold")
        threshold = None if threshold in (None, "") else float(threshold)
        score_type = str(row.get("score_type") or "full").strip().casefold()
        if score_type not in {"full", "domain"}:
            raise ValueError(f"Unsupported HMM score_type {score_type!r} for profile {profile_id!r}")
        requests_by_source.setdefault(source, {})[profile_id] = (canonical, threshold, score_type)

    selected_records: list[tuple[str, str, str, float | None, str]] = []
    for source, requested in requests_by_source.items():
        opener = gzip.open if source.suffix.lower() == ".gz" else open
        record: list[str] = []
        found_for_source: set[str] = set()
        with opener(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                record.append(line)
                if line.strip() != "//":
                    continue
                name = next((value.split(maxsplit=1)[1].strip() for value in record if value.startswith("NAME")), "")
                accession = next((value.split(maxsplit=1)[1].strip().split(".", 1)[0] for value in record if value.startswith("ACC")), "")
                matched_key = next((key for key in requested if key in {name, accession}), None)
                if matched_key is not None:
                    canonical, threshold, score_type = requested[matched_key]
                    selected_records.append(
                        (canonical, matched_key, "".join(record).lstrip(), threshold, score_type)
                    )
                    found_for_source.add(matched_key)
                record = []
        missing_profiles = set(requested) - found_for_source
        if missing_profiles:
            raise ValueError(f"Profiles not found in {source}: {sorted(missing_profiles)}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    cutoff_rows: list[dict[str, Any]] = []
    used_names: set[str] = set()
    with output_path.open("w", encoding="utf-8") as output:
        for canonical, source_profile, record, threshold, score_type in selected_records:
            safe_profile = "".join(character if character.isalnum() or character in "_.-" else "_" for character in source_profile)
            output_name = f"{canonical}__{safe_profile}"
            if output_name in used_names:
                raise ValueError(f"Duplicate selected HMM profile: {output_name}")
            used_names.add(output_name)
            rewritten = []
            for line in record.splitlines(keepends=True):
                rewritten_line = f"NAME  {output_name}" if line.startswith("NAME") else line.rstrip()
                rewritten.append(rewritten_line + "\n")
            output.write("".join(rewritten))
            counts[canonical] = counts.get(canonical, 0) + 1
            if threshold is not None:
                cutoff_rows.append(
                    {"profile_id": output_name, "score_threshold": threshold, "score_type": score_type}
                )

    if cutoff_output is not None:
        cutoff_path = Path(cutoff_output)
        cutoff_path.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(
            cutoff_rows,
            schema={"profile_id": pl.Utf8, "score_threshold": pl.Float64, "score_type": pl.Utf8},
        ).write_csv(cutoff_path)
    return counts


def _confidence(score: float, credible: bool = False) -> str:
    if score <= 0:
        return "none"
    if not credible:
        return "partial"
    if score < 0.7:
        return "low"
    if score < 0.9:
        return "medium"
    return "high"


class MarkerClassifier:
    """Evaluate pathway marker panels and aggregate their COD-group evidence."""

    def __init__(self, catalog: MarkerCatalog | None = None):
        self.catalog = catalog or MarkerCatalog.from_json()

    @staticmethod
    def _require_columns(frame: pl.DataFrame, columns: Iterable[str], table_name: str) -> None:
        missing = set(columns) - set(frame.columns)
        if missing:
            raise ValueError(f"{table_name} is missing required columns: {sorted(missing)}")

    def classify_hits(
        self,
        hits: pl.DataFrame,
        *,
        genome_column: str = "genome_id",
        marker_column: str = "marker_id",
        min_bits: float | None = None,
        max_evalue: float | None = None,
        min_identity: float | None = None,
        min_coverage: float | None = None,
    ) -> pl.DataFrame:
        """Classify genome marker hits into every current e-ADM COD group.

        Thresholds are applied only when explicitly requested.  This lets the
        same API consume pre-filtered annotation presence tables or rawer hits.
        """
        self._require_columns(hits, (genome_column, marker_column), "marker hits")
        genomes = (
            hits.select(pl.col(genome_column).cast(pl.Utf8).alias("genome_id"))
            .drop_nulls()
            .unique(maintain_order=True)["genome_id"]
            .to_list()
        )
        filtered = hits
        threshold_specs = (
            (min_bits, "bits", ">="),
            (max_evalue, "evalue", "<="),
            (min_identity, "identity", ">="),
            (min_coverage, "coverage", ">="),
        )
        for threshold, column, operator in threshold_specs:
            if threshold is None:
                continue
            self._require_columns(filtered, (column,), "marker hits")
            expression = pl.col(column).cast(pl.Float64)
            filtered = filtered.filter(expression >= threshold if operator == ">=" else expression <= threshold)

        presence: dict[str, set[str]] = {genome: set() for genome in genomes}
        for genome, raw_marker in filtered.select(genome_column, marker_column).iter_rows():
            if genome is None or raw_marker is None:
                continue
            canonical = self.catalog.alias_lookup.get(str(raw_marker).strip().casefold())
            if canonical is not None:
                presence.setdefault(str(genome), set()).add(canonical)
        return self.classify_presence(presence)

    def classify_presence(self, genome_markers: Mapping[str, Iterable[str]]) -> pl.DataFrame:
        """Classify a mapping of genome id to canonical marker ids or aliases."""
        rows: list[dict[str, Any]] = []
        for genome_id, raw_markers in genome_markers.items():
            present = {
                canonical
                for marker in raw_markers
                if (canonical := self.catalog.alias_lookup.get(str(marker).strip().casefold())) is not None
            }
            for cod_group in self.catalog.cod_groups:
                group = self.catalog.groups.get(cod_group)
                if group is None:
                    rows.append(
                        self._result_row(
                            str(genome_id), cod_group, False, "", 0.0, False, (),
                            ("no configured rule",), 0, 0,
                        )
                    )
                    continue

                panel_results = [self._evaluate_panel(panel, present) for panel in group.panels]
                best = max(panel_results, key=lambda result: (result[0], result[1], len(result[3]), result[4]))
                score, credible, panel_id, matched, total, missing, weighted, requirements, marker_coverage = best
                rows.append(
                    self._result_row(
                        str(genome_id), cod_group, True, panel_id, score, credible, matched, missing,
                        len(matched), total, weighted, requirements, marker_coverage,
                    )
                )

        schema = {
            "genome_id": pl.Utf8,
            "cod_group": pl.Utf8,
            "score": pl.Float64,
            "credible": pl.Boolean,
            "configured": pl.Boolean,
            "panel_id": pl.Utf8,
            "confidence": pl.Utf8,
            "matched_markers": pl.Utf8,
            "missing_requirements": pl.Utf8,
            "markers_detected": pl.Int64,
            "markers_in_panel": pl.Int64,
            "weighted_completeness": pl.Float64,
            "requirement_coverage": pl.Float64,
            "marker_coverage": pl.Float64,
        }
        return pl.DataFrame(rows, schema=schema)

    def profile_from_counts(
        self,
        marker_counts: Mapping[str, int | float],
        *,
        normalize: bool = True,
        count_weighted: bool = False,
    ) -> dict[str, float]:
        """Convert marker-family counts into a COD-group potential profile.

        Pathway potential is continuous even when the full credibility rule is
        not satisfied. Genome annotations should use the default presence-based
        score. Read-level workflows may request ``count_weighted`` to scale the
        potential by the median observed marker-family count.
        """
        canonical_counts: dict[str, float] = {}
        for raw_marker, raw_count in marker_counts.items():
            canonical = self.catalog.alias_lookup.get(str(raw_marker).strip().casefold())
            count = float(raw_count)
            if canonical is not None and count > 0:
                canonical_counts[canonical] = canonical_counts.get(canonical, 0.0) + count
        classified = self.classify_presence({"sample": canonical_counts})
        profile = {group: 0.0 for group in self.catalog.cod_groups}
        for row in classified.filter(pl.col("score") > 0).to_dicts():
            matched = [marker for marker in row["matched_markers"].split(";") if marker]
            signals = sorted(canonical_counts[marker] for marker in matched if marker in canonical_counts)
            magnitude = 1.0
            if count_weighted:
                if not signals:
                    continue
                midpoint = len(signals) // 2
                magnitude = signals[midpoint] if len(signals) % 2 else (signals[midpoint - 1] + signals[midpoint]) / 2
            profile[row["cod_group"]] = float(row["score"]) * magnitude
        total = sum(profile.values())
        if normalize and total > 0:
            profile = {group: value / total for group, value in profile.items()}
        return profile

    def community_profile_from_hits(
        self,
        hits: pl.DataFrame,
        genome_abundances: Mapping[str, int | float],
        *,
        normalize: bool = False,
        copy_cap: float = 1.0,
    ) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, float]]:
        """Aggregate genome marker evidence before calculating COD potential.

        Assumptions
        -----------
        * Genome abundances are fractions of the complete sample when their
          positive total is at most one. Count-like inputs whose total exceeds
          one are converted to fractions.
        * A marker family contributes at most once per genome by default. This
          avoids treating assembly fragmentation, paralogs, or uncertain gene
          copy number as proportional metabolic activity. ``copy_cap`` keeps
          this policy explicit and extensible.
        * Accepted marker hits are pooled across the community before pathway
          rules are evaluated. Panel requirements and credibility thresholds
          remain diagnostic; they do not zero community potential.
        * Marker weights are normalized within each panel. Alternative panels
          represent alternative biochemical routes, so the best-supported
          panel defines a COD group's raw potential. This prevents groups with
          more catalogued alternatives from receiving a larger prior weight.
        * The returned potential describes relative gene-encoded capacity, not
          expression, flux, biomass yield, or measured activity.
        """
        self._require_columns(hits, ("genome_id", "gene_id", "marker_id"), "marker hits")
        if copy_cap <= 0:
            raise ValueError("copy_cap must be positive")

        positive_abundances = {
            str(genome): float(abundance)
            for genome, abundance in genome_abundances.items()
            if float(abundance) > 0
        }
        abundance_total = sum(positive_abundances.values())
        abundance_scale = 1.0 if abundance_total <= 1.0 + 1e-9 else abundance_total
        scaled_abundances = {
            genome: abundance / abundance_scale
            for genome, abundance in positive_abundances.items()
        }

        gene_ids: dict[tuple[str, str], set[str]] = {}
        for genome, gene, raw_marker in hits.select("genome_id", "gene_id", "marker_id").iter_rows():
            genome = str(genome)
            if genome not in scaled_abundances or gene is None or raw_marker is None:
                continue
            canonical = self.catalog.alias_lookup.get(str(raw_marker).strip().casefold())
            if canonical is None:
                continue
            gene_ids.setdefault((genome, canonical), set()).add(str(gene))

        marker_weighted_abundance = {marker: 0.0 for marker in self.catalog.markers}
        marker_genomes = {marker: 0 for marker in self.catalog.markers}
        marker_gene_copies = {marker: 0 for marker in self.catalog.markers}
        for (genome, marker), genes in gene_ids.items():
            copies = len(genes)
            marker_genomes[marker] += 1
            marker_gene_copies[marker] += copies
            marker_weighted_abundance[marker] += scaled_abundances[genome] * min(float(copies), copy_cap)

        marker_rows = []
        for marker, definition in self.catalog.markers.items():
            marker_rows.append(
                {
                    "marker_id": marker,
                    "marker_name": str(definition.get("name", marker)),
                    "genomes_detected": marker_genomes[marker],
                    "gene_copies_detected": marker_gene_copies[marker],
                    "weighted_abundance": marker_weighted_abundance[marker],
                    "copy_cap_per_genome": float(copy_cap),
                }
            )
        marker_frame = pl.DataFrame(
            marker_rows,
            schema={
                "marker_id": pl.Utf8,
                "marker_name": pl.Utf8,
                "genomes_detected": pl.Int64,
                "gene_copies_detected": pl.Int64,
                "weighted_abundance": pl.Float64,
                "copy_cap_per_genome": pl.Float64,
            },
        )

        community_present = {
            marker for marker, abundance in marker_weighted_abundance.items() if abundance > 0
        }
        panel_rows: list[dict[str, Any]] = []
        profile = {group: 0.0 for group in self.catalog.cod_groups}
        for cod_group in self.catalog.cod_groups:
            group = self.catalog.groups.get(cod_group)
            if group is None:
                continue
            group_rows: list[dict[str, Any]] = []
            for panel in group.panels:
                total_weight = sum(panel.marker_weights.values())
                contributions = {
                    marker: marker_weighted_abundance[marker] * weight / total_weight
                    for marker, weight in panel.marker_weights.items()
                    if marker_weighted_abundance[marker] > 0
                }
                potential = sum(contributions.values())
                _, credible, _, matched, total, missing, weighted, requirements, coverage = (
                    self._evaluate_panel(panel, community_present)
                )
                group_rows.append(
                    {
                        "cod_group": cod_group,
                        "panel_id": panel.id,
                        "potential": potential,
                        "selected": False,
                        "credible": credible,
                        "matched_markers": ";".join(matched),
                        "missing_requirements": ";".join(missing),
                        "markers_detected": len(matched),
                        "markers_in_panel": total,
                        "weighted_completeness": weighted,
                        "requirement_coverage": requirements,
                        "marker_coverage": coverage,
                        "marker_contributions": ";".join(
                            f"{marker}:{value:.12g}" for marker, value in sorted(contributions.items())
                        ),
                    }
                )
            selected = max(group_rows, key=lambda row: (row["potential"], row["credible"], row["markers_detected"]))
            selected["selected"] = True
            profile[cod_group] = float(selected["potential"])
            panel_rows.extend(group_rows)

        panel_frame = pl.DataFrame(
            panel_rows,
            schema={
                "cod_group": pl.Utf8,
                "panel_id": pl.Utf8,
                "potential": pl.Float64,
                "selected": pl.Boolean,
                "credible": pl.Boolean,
                "matched_markers": pl.Utf8,
                "missing_requirements": pl.Utf8,
                "markers_detected": pl.Int64,
                "markers_in_panel": pl.Int64,
                "weighted_completeness": pl.Float64,
                "requirement_coverage": pl.Float64,
                "marker_coverage": pl.Float64,
                "marker_contributions": pl.Utf8,
            },
        )
        total = sum(profile.values())
        if normalize and total > 0:
            profile = {group: value / total for group, value in profile.items()}
        return marker_frame, panel_frame, profile

    @staticmethod
    def _evaluate_panel(
        panel: MarkerPanel,
        present: set[str],
    ) -> tuple[float, bool, str, tuple[str, ...], int, tuple[str, ...], float, float, float]:
        matched = tuple(sorted(set(panel.marker_weights) & present))
        missing: list[str] = [marker for marker in panel.required_all if marker not in present]
        requirements_met = sum(marker in present for marker in panel.required_all)
        for clause in panel.required_any:
            if not set(clause) & present:
                missing.append("one_of(" + "|".join(clause) + ")")
            else:
                requirements_met += 1
        weighted_score = sum(panel.marker_weights[marker] for marker in matched) / sum(panel.marker_weights.values())
        requirement_total = len(panel.required_all) + len(panel.required_any)
        requirement_coverage = requirements_met / requirement_total if requirement_total else 1.0
        marker_coverage = min(1.0, len(matched) / panel.minimum_markers)
        requirements_pass = not missing and len(matched) >= panel.minimum_markers
        credible = requirements_pass and weighted_score >= panel.minimum_score
        if not matched or (panel.strict and not credible):
            score = 0.0
        else:
            score = weighted_score * (0.35 + 0.65 * requirement_coverage) * (0.35 + 0.65 * marker_coverage)
        return (
            score, credible, panel.id, matched, len(panel.marker_weights), tuple(missing),
            weighted_score, requirement_coverage, marker_coverage,
        )

    @staticmethod
    def _result_row(
        genome_id: str,
        cod_group: str,
        configured: bool,
        panel_id: str,
        score: float,
        credible: bool,
        matched: Iterable[str],
        missing: Iterable[str],
        detected: int,
        total: int,
        weighted_completeness: float = 0.0,
        requirement_coverage: float = 0.0,
        marker_coverage: float = 0.0,
    ) -> dict[str, Any]:
        return {
            "genome_id": genome_id,
            "cod_group": cod_group,
            "score": score,
            "credible": credible,
            "configured": configured,
            "panel_id": panel_id,
            "confidence": _confidence(score, credible),
            "matched_markers": ";".join(matched),
            "missing_requirements": ";".join(missing),
            "markers_detected": detected,
            "markers_in_panel": total,
            "weighted_completeness": weighted_completeness,
            "requirement_coverage": requirement_coverage,
            "marker_coverage": marker_coverage,
        }

    def aggregate_samples(
        self,
        genome_scores: pl.DataFrame,
        abundances: pl.DataFrame,
        *,
        sample_column: str = "sample",
        genome_column: str = "genome_id",
        abundance_column: str = "abundance",
        normalize: bool = True,
    ) -> tuple[pl.DataFrame, pl.DataFrame]:
        """Abundance-weight genome evidence into sample-by-COD-group profiles."""
        self._require_columns(genome_scores, ("genome_id", "cod_group", "score"), "genome scores")
        self._require_columns(abundances, (sample_column, genome_column, abundance_column), "abundances")

        abundance = (
            abundances.select(
                pl.col(sample_column).cast(pl.Utf8).alias("sample"),
                pl.col(genome_column).cast(pl.Utf8).alias("genome_id"),
                pl.col(abundance_column).cast(pl.Float64).alias("abundance"),
            )
            .drop_nulls()
            .filter(pl.col("abundance") >= 0)
            .group_by("sample", "genome_id")
            .agg(pl.col("abundance").sum())
        )
        if abundance.height == 0:
            raise ValueError("abundances contains no non-negative numeric values")
        totals = abundance.group_by("sample").agg(pl.col("abundance").sum().alias("total_abundance"))
        if totals.filter(pl.col("total_abundance") <= 0).height:
            raise ValueError("Every sample must have positive total genome abundance")
        abundance = abundance.join(totals, on="sample").with_columns(
            pl.when(pl.col("total_abundance") > 1.0 + 1e-9)
            .then(pl.col("abundance") / pl.col("total_abundance"))
            .otherwise(pl.col("abundance"))
            .alias("genome_weight")
        )

        scores = genome_scores.select(
            pl.col("genome_id").cast(pl.Utf8),
            pl.col("cod_group").cast(pl.Utf8),
            pl.col("score").cast(pl.Float64),
        ).group_by("genome_id", "cod_group").agg(pl.col("score").max())
        weighted = (
            abundance.join(scores, on="genome_id", how="left")
            .with_columns(pl.col("score").fill_null(0.0))
            .with_columns((pl.col("genome_weight") * pl.col("score")).alias("weighted_score"))
            .group_by("sample", "cod_group")
            .agg(pl.col("weighted_score").sum().alias("abundance"))
        )
        samples = abundance.select("sample").unique(maintain_order=True)
        cod_groups = pl.DataFrame({"cod_group": list(self.catalog.cod_groups)})
        profile = samples.join(cod_groups, how="cross").join(weighted, on=["sample", "cod_group"], how="left")
        profile = profile.with_columns(pl.col("abundance").fill_null(0.0))
        if normalize:
            profile = profile.with_columns(pl.col("abundance").sum().over("sample").alias("_cod_total")).with_columns(
                pl.when(pl.col("_cod_total") > 0)
                .then(pl.col("abundance") / pl.col("_cod_total"))
                .otherwise(0.0)
                .alias("abundance")
            ).drop("_cod_total")

        genome_max = scores.group_by("genome_id").agg(pl.col("score").max().alias("max_score"))
        qc = (
            abundance.join(genome_max, on="genome_id", how="left")
            .with_columns(pl.col("max_score").fill_null(0.0))
            .group_by("sample")
            .agg(
                pl.col("total_abundance").first(),
                pl.when(pl.col("max_score") > 0).then(pl.col("abundance")).otherwise(0.0).sum().alias("classified_abundance"),
            )
            .with_columns((pl.col("classified_abundance") / pl.col("total_abundance")).alias("classified_fraction"))
        )
        detected = profile.filter(pl.col("abundance") > 0).group_by("sample").agg(pl.len().alias("groups_detected"))
        qc = qc.join(detected, on="sample", how="left").with_columns(pl.col("groups_detected").fill_null(0))
        return profile.sort("sample", "cod_group"), qc.sort("sample")
