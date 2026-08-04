# Gene-Marker Biochemistry

This page documents the **biochemistry** behind ADToolbox's direct gene-marker
allocator: which anaerobic-digestion function each e-ADM COD group represents,
the metabolic pathway that defines it, the diagnostic gene families that mark
that pathway, and the rules that turn marker presence into a COD-group score.

For the *workflow* — how to run the allocator, build marker databases, choose
the MMseqs/HMMER backend, and read the output files — see
[Direct gene-marker allocation](Metagenomics_Pipeline.md#direct-gene-marker-allocation).
This page is the companion *reference* for the rules themselves.

- **Catalog file:** `adtoolbox/pkg_data/gene_marker_catalog.json`
- **Catalog version:** `0.3.0` (`catalog_version` field; also stamped into cache filenames)
- **COD groups covered:** all 15 microbial e-ADM groups, each with ≥1 curated panel
- **Every marker carries an `evidence` field** — a KEGG/UniProt/PMC URL for the
  gene family — so each assignment is traceable to a source.

> **What a score means.** Marker DNA establishes *encoded functional potential*,
> not expression, flux, growth rate, reaction *direction*, or measured activity.
> A high score says "this community encodes the machinery for this conversion,"
> not "this conversion is running." Keep this in mind especially for the
> reversible reactions (chain elongation vs. syntrophic β-oxidation).

---

## Files that define the biochemistry

> **Is it a TOML file? No.** The biochemistry lives in a single **JSON** rule
> catalog. TOML *is* used elsewhere in ADToolbox — the pipeline's
> [execution profiles](Metagenomics_Pipeline.md#execution-profiles) are TOML
> (see that section for a runtime example and the full reference profile) —
> but those configure Slurm and step settings, not the biochemistry.

The biochemistry is split into **two deliberately separate layers** plus the
interpreter code:

1. **The rules** — which gene families imply which COD group. One JSON file.
2. **The detection assets** — the sequences/profiles used to *find* those gene
   families in a genome or read set (FASTA for MMseqs, HMMs for HMMER).

They are kept apart on purpose: changing a scoring rule must not silently change
sequence homology, and adding reference diversity must not change the biochemical
interpretation.

| Layer | File | Format | Role |
| --- | --- | --- | --- |
| **Rules** | `adtoolbox/pkg_data/gene_marker_catalog.json` | JSON | **The biochemistry.** 96 markers (with aliases + `evidence` URLs), 15 COD groups, 25 panels, and every threshold. This is the file you edit to change the biochemistry. |
| Interpreter | `adtoolbox/markers.py` | Python | Loads and validates the catalog, scores panels, parses MMseqs/HMMER hits. Contains *no* gene→group knowledge — that all lives in the JSON. |
| Detection (HMMER) | `…/marker_databases/v0.2.0/Marker_Profiles.hmm` | HMMER3 | Profile HMMs (selected KOfam + dbCAN families). |
| Detection (HMMER) | `…/marker_databases/v0.2.0/Marker_Profile_Cutoffs.csv` | CSV | Per-profile adaptive score thresholds and full/best-domain score type. |
| Detection (MMseqs) | `…/marker_databases/v0.2.0/Marker_Protein_DB.fasta` | FASTA | One `hmmemit` **consensus** sequence per profile; headers are `<accession>\|<marker_id>`. Not isolate proteins. |
| Detection (MMseqs) | `…/marker_databases/v0.2.0/Marker_Protein_DB_mmseqs*` | MMseqs DB | Prebuilt searchable target, reused across samples. |
| Provenance | `…/marker_databases/v0.2.0/profile_manifest.tsv` | TSV | Auditable marker → source-profile mapping. |
| Provenance | `…/marker_databases/v0.2.0/SOURCES.json` | JSON | Pinned upstream KOfam/dbCAN releases, URLs, and license notes. |
| Provenance | `…/marker_databases/v0.2.0/SHA256SUMS` | text | Integrity checks for the shipped assets. |
| Provenance | `…/marker_databases/v0.2.0/README.md` | Markdown | Human description of the bundled assets. |
| Rebuild | `…/marker_databases/v0.2.0/rebuild.py` | Python | Regenerates every detection asset from the catalog + upstream sources. |

All of these are auto-resolved: `MarkerCatalog.from_json()` defaults to the
bundled catalog, and the detection assets are selected automatically when
`--protein-db` / `--marker-hmm-db` are omitted.

> **Version gotcha.** The directory is named `v0.2.0` but its **contents are
> catalog `0.3.0`** (see the `README.md` title and `SOURCES.json → asset_version`).
> The stable directory name is retained on purpose so existing installations keep
> resolving the bundled paths; the *catalog* version (`0.3.0`) is the one stamped
> into COD-cache filenames.

---

## How a panel is scored

Each COD group owns one or more **panels**. A panel is one biochemical route to
the group's function (e.g. `X_ch` has separate cellulose, xylan, starch, and
pectin panels). A genome is scored against every panel; the **best-supported
panel** supplies the group's score. Panels are *alternatives*, never additive —
a group with more catalogued routes does not get a larger prior.

A panel declares:

| Field | Meaning |
| --- | --- |
| `marker_weights` | The gene families in the panel and their diagnostic weight (higher = more specific to this function). |
| `required_all` | Markers that **must all** be present for the pathway to be credible. |
| `required_any` | A list of clauses; **each clause** must contribute ≥1 present marker (an "AND of ORs"). |
| `minimum_markers` | Minimum number of distinct panel markers that must be detected. |
| `minimum_score` | Minimum *weighted completeness* (0–1) for the panel to be called credible. |
| `strict` | If `true`, a non-credible panel scores **0** instead of a partial score. Reserved for specialized functions where a false positive would badly mislead the model. |

From these, the classifier computes:

- **weighted completeness** = (Σ weights of matched markers) / (Σ all panel weights)
- **requirement coverage** = fraction of `required_all` + `required_any` clauses satisfied
- **marker coverage** = min(1, matched / `minimum_markers`)
- **credible** = all requirements met **and** ≥ `minimum_markers` **and** weighted completeness ≥ `minimum_score`
- **score** = `weighted_completeness × (0.35 + 0.65·requirement_coverage) × (0.35 + 0.65·marker_coverage)`
  — unless the panel is `strict` and not credible, in which case the score is 0.

The per-genome `confidence` tier is `none` / `partial` / `low` / `medium` /
`high`, derived from the score and the credibility flag.

> **Per-genome vs. community.** The credibility rule above is a strict
> *per-genome* diagnostic. The **final** COD profile is built from the
> abundance-weighted *community marker pool*, which is deliberately more
> permissive: accepted marker evidence from the community is not zeroed just
> because no single genome encoded the complete textbook pathway. See
> [Community marker-pool COD potential and credibility](Metagenomics_Pipeline.md#community-marker-pool-cod-potential-and-credibility).

---

## COD groups by digestion stage

### 1 · Hydrolysis of particulates

#### `X_ch` — carbohydrate hydrolysis

Extracellular depolymerization of particulate polysaccharides to sugars. Four
alternative panels cover the major substrates; each pairs catalytic CAZyme
families with binding, transport, or downstream-utilization evidence so that a
lone glycosidase cannot classify a genome.

| Panel | Pathway | Key markers |
| --- | --- | --- |
| `cellulose_utilization` | Cellulose → cellobiose → glucose | `GH5`/`GH9`/`GH48` endo/processive cellulases, `CBM3` cellulose-binding module, `bgl` β-glucosidase (required) |
| `xylan_utilization` | Xylan → xylose → catabolism | `GH10`/`GH11` endoxylanases, `GH43` debranching, `CBM6`, plus `xylA`/`xylB` (xylose isomerase + xylulokinase, required) |
| `starch_utilization` | Starch → maltodextrin → glucose | `GH13` amylase/pullulanase (required), `agl` α-glucosidase, SusC/SusD glycan uptake |
| `pectin_utilization` | Pectin de-esterification + cleavage | `PL1` pectate lyase (required), `CE8` pectin methylesterase, `GH43`, SusC/SusD |

#### `X_pr` — protein hydrolysis

Secreted/cell-surface proteolysis coupled to oligopeptide acquisition. Two
panels: `proteolysis_and_peptide_import` (a secreted endopeptidase **and** a
peptide transporter **and** an aminopeptidase) and the broader
`peptide_scavenging` panel.

Markers: `extracellular_protease` (subtilisin/aprE/prtP), `endopeptidase`,
aminopeptidases (`pepN`, `pepO`, `pepA`, `pepT`, `pepQ`, `pepD`), and the
oligo/di-peptide transporters `oppA/oppB`, `dppA/dppB/dppC`, `dtpT`.

#### `X_li` — lipid hydrolysis

Extracellular ester hydrolysis with fatty-acid activation. Markers: `lipase`
(triacylglycerol lipase), `esterase`, `phospholipase`, `lysophospholipase`,
`acylglycerol_lipase`, and `fadD` (long-chain fatty-acyl-CoA ligase) linking
hydrolysis to β-oxidation.

### 2 · Acidogenesis (monomer fermentation)

#### `X_su` — sugar fermentation

Uptake and fermentation of monosaccharides to pyruvate and mixed acids. Each
panel requires a **transporter**, **central-pathway coverage**, and a
**fermentative outlet** — generic glycolysis alone does not qualify.

| Panel | Central pathway | Required backbone |
| --- | --- | --- |
| `hexose_fermentation` | Embden–Meyerhof–Parnas (EMP) | `pgi`, `pfkA`, `gapA`, `pyk` + transporter (`ptsG`/`glcP`) + outlet (`pfor`/`adhE`/`lctD`) |
| `pentose_fermentation` | Pentose → EMP via isomerase | `xylA`, `xylB`, `gapA`, `pyk` + `xylF` + outlet |
| `entner_doudoroff` | Entner–Doudoroff (ED) | `edd`, `eda` + transporter + `gapA`/`pyk` + outlet |

#### `X_aa` — amino-acid fermentation

Diagnostic Stickland acceptor pathways plus a general catabolism panel.

- `glycine_stickland`: glycine reductase `grdA`/`grdB` (required) + peptide
  import + `pfor`.
- `proline_stickland`: D-proline reductase `prdA`/`prdB` (required) + peptide
  import + `pfor`.
- `general_amino_acid_catabolism`: `gdh` (glutamate dehydrogenase), `ilvE`
  (branched-chain aminotransferase), `tdcB` (threonine dehydratase), `pfor`,
  and a peptide transporter.

The Stickland reductases are used as *diagnostic* acceptor-side markers because
they are specific to anaerobic amino-acid fermentation, unlike the widely
distributed transaminases.

#### `X_fa` — long-chain fatty-acid β-oxidation

Activation and β-oxidation of LCFAs to acetyl-CoA. Single panel requiring the
full cycle: `fadD` (activation), `fadE`/`acd` (acyl-CoA dehydrogenase), `fadB`
(hydratase/dehydrogenase), `fadA` (thiolase), plus electron-transfer
flavoprotein `etfA`/`etfB`.

### 3 · Secondary fermentation & chain elongation

#### `X_et` — ethanol oxidation

Ethanol → acetaldehyde → acetyl-CoA → acetate, **without** requiring
chain-elongation capacity. Markers: `adhE` (bifunctional), `adh`, `aldh`, and
the acetate node `acs`/`pta`/`ackA`.

#### `X_lac` — lactate oxidation

Lactate → pyruvate via the electron-confurcating lactate dehydrogenase
`lctD`(+`lctB`/`lctC`) and `pfor`, with `rnfB` for energy conservation. Required:
`lctD` + `pfor`.

#### `X_ac_et` / `X_ac_lac` — acetate-consuming chain initiation

Chain elongation that draws acetyl-CoA from **acetate** using an
electron-donor-derived reduced pool (ethanol for `X_ac_et`, lactate for
`X_ac_lac`). Both panels are `strict`. Shared reverse-β-oxidation core: `thl`
(thiolase), `hbd`, `crt`, `bcd`, `etfA`/`etfB`. Donor-specific: ethanol
(`adhE`/`adh`/`aldh`) or lactate (`lctD`/`lctB`/`lctC`/`pfor`). Acetate node:
`acs`/`ackA`/`pta`.

#### `X_chain_et` / `X_chain_lac` — reverse β-oxidation chain elongation

Elongation of VFAs by C2 units to medium-chain acids (butyrate → caproate),
driven by ethanol or lactate. Both `strict`. Core `bcd`/`etfA`/`etfB` plus the
**terminal CoA transferase `cat`/`but`** (required) — the release step that
distinguishes elongation from the oxidative direction (see below).

### 4 · Syntrophic VFA oxidation

#### `X_VFA_deg` — syntrophic butyrate & propionate oxidation

Oxidation of VFAs to acetate/H₂/CO₂ under low-H₂ syntrophy. Two panels:

- `syntrophic_butyrate_oxidation`: the β-oxidation core (`bcd`/`etf` + `thl`/
  `hbd`/`crt`) **plus** interspecies electron-transfer evidence
  `hydA`/`fdhA`/`rnfB`.
- `syntrophic_propionate_mmc`: the **methylmalonyl-CoA (MMC) pathway** —
  `mmcE`/`mmcF` (mutase), `mmcG` (epimerase), `sucC`/`sucD` (succinyl-CoA
  synthetase), `fumB`, `mdh` — plus `hydA`/`fdhA`.

Not `strict`, because these enzymes double as the biosynthetic/oxidative
direction; the hydrogen/formate-transfer markers raise confidence but cannot by
themselves prove the *in situ* direction.

### 5 · Methanogenesis

#### `X_Me_ac` — acetoclastic methanogenesis

Acetate → CH₄ + CO₂. Two `strict` panels differing in the acetate-activation
route, both requiring `mcrA` (methyl-CoM reductase) and a CODH/ACS subunit
(`cdhC`/`cdhD`/`cdhE`):

- `acetoclastic_acs`: AMP-forming `acs` activation (*Methanothrix*/*Methanosaeta*
  type — the dominant acetoclastic genus in digesters).
- `acetoclastic_ack_pta`: `ackA` + `pta` activation (*Methanosarcina* type).

#### `X_Me_CO2` — hydrogenotrophic methanogenesis

CO₂ + H₂ → CH₄ via the H₄MPT C1-carrier pathway. Single `strict` panel:
`mcrA` (required) plus the CO₂-reduction branch `fwdA` (formyl-MF
dehydrogenase), `ftr`, `mch`, `mtd`/`mer`, with `mtrA` (methyltransferase).

---

## Direction, ambiguity, and known limitations

- **Reverse β-oxidation is reversible.** The same `thl`/`hbd`/`crt`/`bcd`/`etf`
  gene set runs *forward* for syntrophic VFA oxidation (`X_VFA_deg`) and
  *reverse* for chain elongation (`X_chain_*`). The catalog leans elongation vs.
  oxidation on the terminal transferase `cat`/`but` (elongation) vs.
  hydrogen/formate transfer `hydA`/`fdhA` (syntrophy), but **DNA cannot fix the
  active direction** — feed chemistry or expression data is needed to resolve it.
- **`strict` groups are all-or-nothing per genome.** The specialized groups
  (`X_Me_ac`, `X_Me_CO2`, `X_ac_*`, `X_chain_*`) score 0 for a genome unless the
  panel is fully credible. This is intentional (a false methanogen is worse than
  a missed one), and the community marker pool still credits partial evidence.
- **Methylotrophic methanogenesis is not yet catalogued.** The two methanogenesis
  groups cover the acetoclastic and hydrogenotrophic routes; methanol/methylamine
  (`mta`/`mtb`/`mtt`) methanogens are not represented and will not score in
  either group.
- **Markers are specific orthologs.** Organisms using non-orthologous isoenzymes
  for a step, or dividing a pathway across community members (cross-feeding), may
  be under-scored per genome — another reason the community pool is the reporting
  unit rather than the single-genome credibility call.

---

## Editing or extending the catalog

The rules are **data, not code** — you never touch Python to change the
biochemistry. The catalog has two top-level sections you edit:

- `markers` — the dictionary of gene families. Each entry has a canonical id, a
  list of `aliases` (KO numbers, gene symbols, CAZy families — all matched
  case-insensitively), a human `name`, and an `evidence` URL.
- `groups` — one entry per COD group, each holding a list of `panels`
  (alternative biochemical routes). A panel is the unit that gets scored.

The general loop is: **copy → edit → validate → (optionally add detection
sequences) → run with `--catalog`.**

```bash
cp adtoolbox/pkg_data/gene_marker_catalog.json my_catalog.json
# ...edit my_catalog.json...
adtoolbox metagenomics validate-marker-catalog --catalog my_catalog.json
adtoolbox metagenomics process ... --catalog my_catalog.json
```

`validate-marker-catalog` is strict — it rejects unknown COD groups, empty
`required_any` clauses, requirement markers that have no weight, `minimum_markers`
outside `1…len(panel)`, `minimum_score` outside `0…1`, and duplicate aliases —
so most mistakes are caught before a run. The catalog version is stamped into
COD-cache filenames, so bumping `catalog_version` guarantees stale profiles are
not silently reused.

### Example 1 — Retune a threshold (loosen or tighten a panel)

The simplest edit changes only the numbers/requirements on an existing panel.
Say you want the hydrogenotrophic-methanogenesis panel to accept a genome on
`mcrA` plus **any two** C1-pathway genes, rather than demanding the full pathway.
Edit its `required_all`/`required_any`, lower `minimum_markers`, and (optionally)
drop `strict`:

```jsonc
// groups → X_Me_CO2 → panels[0]
{
  "id": "hydrogenotrophic_core",
  "marker_weights": {"mcrA": 2, "mcrB": 1, "mtrA": 1,
                     "fwdA": 1, "ftr": 1, "mch": 1, "mtd": 1, "mer": 1},
  "required_all": ["mcrA"],
  "required_any": [["fwdA", "ftr", "mch", "mtd", "mer"]],  // was 4 separate clauses
  "minimum_markers": 3,                                      // was 6
  "minimum_score": 0.45,                                     // was 0.65
  "strict": false                                            // was true
}
```

Nothing else changes — no new markers, no database edit — because every marker
referenced already exists and is already in the detection assets.

### Example 2 — Add an alternative panel to an existing group

Panels are alternatives, so you can add a second route to a group without
touching the first. To recognize *methylotrophic-leaning* acetate handling via a
different activation, append a panel to `X_Me_ac.panels` (all its markers must
already exist in `markers`):

```jsonc
// groups → X_Me_ac → panels  (append)
{
  "id": "acetoclastic_minimal",
  "marker_weights": {"mcrA": 2, "mcrB": 1, "acs": 2},
  "required_all": ["mcrA", "acs"],
  "required_any": [],
  "minimum_markers": 2,
  "minimum_score": 0.6,
  "strict": true
}
```

The best-supported panel wins, so adding a route can only *increase* a group's
score for genomes that match it, never penalize genomes matched by another panel.

### Example 3 — Add a brand-new marker family

Suppose you want to cover **methylotrophic methanogenesis** (a documented gap —
see [limitations](#direction-ambiguity-and-known-limitations)). That needs new
markers *and* a new panel/group. First add the gene families to `markers`:

```jsonc
// markers  (add entries)
"mtaB": {"aliases": ["K14080"], "name": "methanol--corrinoid methyltransferase",
         "evidence": "https://www.genome.jp/entry/K14080"},
"mtbA": {"aliases": ["K14082"], "name": "methylamine--corrinoid methyltransferase",
         "evidence": "https://www.genome.jp/entry/K14082"}
```

Then give them a panel. If you want a new COD group you must also add its name to
the top-level `cod_groups` list (and the e-ADM model must have a matching group);
to stay within the existing 15 groups, attach the panel to `X_Me_CO2` instead:

```jsonc
// groups → X_Me_CO2 → panels  (append)
{
  "id": "methylotrophic",
  "marker_weights": {"mcrA": 2, "mtaB": 2, "mtbA": 2},
  "required_all": ["mcrA"],
  "required_any": [["mtaB", "mtbA"]],
  "minimum_markers": 2,
  "minimum_score": 0.5,
  "strict": true
}
```

Because the score comes from a genome's marker hits, **a new marker does nothing
until the detection database can find it** (Example 4).

### Example 4 — Make a new marker detectable

New markers must be added to whichever backend you run. Both take a manifest, so
you never hand-edit FASTA headers or HMM `NAME` lines.

**MMseqs backend** — collect curated reference proteins and build the target
FASTA (headers are rewritten to `<accession>|<marker_id>`):

```text
# marker_references.tsv
marker_id	source_fasta
mtaB	references/mtaB.faa
mtbA	references/mtbA.faa.gz
```

```bash
adtoolbox database build-marker-protein-db \
  --manifest marker_references.tsv \
  --output My_Marker_Protein_DB.fasta \
  --catalog my_catalog.json

adtoolbox metagenomics validate-marker-protein-db \
  --protein-db My_Marker_Protein_DB.fasta \
  --catalog my_catalog.json          # confirms every catalog marker is covered
```

**HMMER backend** — select profiles from KOfam/dbCAN/Pfam/custom collections and
carry their score thresholds:

```text
# marker_hmms.tsv
marker_id	source_hmm	profile_id	score_threshold	score_type
mtaB	kofam/profiles/K14080.hmm	K14080	120.5	full
mtbA	kofam/profiles/K14082.hmm	K14082	118.2	full
```

```bash
adtoolbox database build-marker-hmm-db \
  --manifest marker_hmms.tsv \
  --output My_Marker_Profiles.hmm \
  --cutoffs-output My_Marker_Profile_Cutoffs.csv \
  --catalog my_catalog.json
```

Then run with your catalog and database together:

```bash
# MMseqs
adtoolbox metagenomics process ... \
  --catalog my_catalog.json \
  --protein-db My_Marker_Protein_DB.fasta

# HMMER
adtoolbox metagenomics process ... \
  --catalog my_catalog.json \
  --marker-backend hmmer \
  --marker-hmm-db My_Marker_Profiles.hmm \
  --marker-hmm-cutoffs My_Marker_Profile_Cutoffs.csv
```

To regenerate the *entire* bundled asset set from the catalog and pinned upstream
sources instead of adding one family, use the shipped rebuild script:

```bash
python adtoolbox/pkg_data/marker_databases/v0.2.0/rebuild.py --help
```

### Editing checklist

- [ ] Every marker in a panel's `marker_weights`, `required_all`, and
  `required_any` exists in the `markers` section.
- [ ] Requirement markers also have a weight (the validator enforces this).
- [ ] `minimum_markers` ≤ number of markers in the panel; `0 ≤ minimum_score ≤ 1`.
- [ ] New COD-group names are added to the top-level `cod_groups` **and** exist in
  the e-ADM model.
- [ ] New markers are added to the detection database (Example 4) and
  `validate-marker-protein-db` reports full coverage.
- [ ] `catalog_version` bumped so caches refresh.
- [ ] Each new marker keeps an `evidence` URL for traceability.

## See also

- [Metagenomics Pipeline — Direct gene-marker allocation](Metagenomics_Pipeline.md#direct-gene-marker-allocation) — the end-to-end workflow and output files.
- [`adtoolbox.markers` API](api-markers.md) — generated reference for `MarkerCatalog`, `MarkerClassifier`, and the hit parsers.
- [ADM Models](ADM_Models.md) — what the COD groups mean inside e-ADM.
