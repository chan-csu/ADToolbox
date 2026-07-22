# Bioconda recipe for ADToolbox

`meta.yaml` here is a standard `noarch: python` Bioconda recipe. Drop it into
`bioconda-recipes/recipes/adtoolbox/meta.yaml` and open a PR as usual.

## Blockers to clear first

These are the only things stopping the recipe from building green:

1. **Publish 1.1.14 to PyPI.** The last published sdist (1.1.13) still does
   `from distutils.log import warn`, which raises `ModuleNotFoundError` on
   Python ≥ 3.12 — so a package built from it fails Bioconda's import test on
   modern interpreters. The fix is in the working tree (`core.py` now uses
   `from warnings import warn`). Publish 1.1.14 from the current tree, then set
   the recipe `sha256`:

   ```bash
   curl -s https://pypi.org/pypi/adtoolbox/json \
     | python -c "import json,sys; d=json.load(sys.stdin); v=d['info']['version']; \
         print(v, [f['digests']['sha256'] for f in d['releases'][v] if f['packagetype']=='sdist'][0])"
   ```

2. **Add a license.** The repo has no `LICENSE` file and no `license` in
   `pyproject.toml`/PyPI/GitHub, so `about.license` is set to `LicenseRef-UNSET`.
   The linter rejects that. Add a `LICENSE` and set `license:` + `license_file:`.

3. **Fix the version string** (recommended). `adtoolbox/__init__.py` hardcodes
   `__version__ = "1.1.0"`, so `adtoolbox --version` prints `1.1.0` even though
   the package is 1.1.14. Derive it from installed metadata instead, or bump it,
   before cutting the release.

## Decision: how heavy should the package be?

`run:` currently pulls the full external-tool stack — `fastp`, `cutadapt`,
`bioconductor-dada2`, `vsearch`, `mmseqs2`, `sra-tools`, `ncbi-datasets-cli`.
That means **every `conda install adtoolbox` drags in R** (via `bioconductor-dada2`)
and the whole bioinformatics stack, even for someone who only runs ADM
simulations.

Two reasonable options:

- **Keep them** — one install gives a fully working pipeline (matches "all
  dependencies figured out").
- **Split** — ship a light `adtoolbox` (Python deps only) and let users add the
  tools via `environment.yml` or a `conda install` of the pipeline tools. Some
  reviewers prefer this for a Python package whose tools are optional
  subprocesses.

Either is fine; it's a packaging-philosophy call. The recipe is written for the
first; delete the tool lines from `run:` (and the tool `--version` lines from
`test:`) for the second.

## Channel notes (verified against the Anaconda API)

- All tools are on **bioconda** except `ncbi-datasets-cli`, which is on
  **conda-forge**.
- `openbox` and `dash-escher` are on **neither** — they stay pip-only extras.
- `dash`, `dash-bootstrap-components`, `pygad`, `pytorch` are on conda-forge if
  you want the dashboard / optimizer extras added to `run:`.
- Don't reintroduce the `pyproject.toml` upper caps (`numpy<2`, `polars<0.21`,
  `rich<13`). They're stale — the suite passes on numpy 2.3, polars 1.41, rich
  14.2 — and copying them would force conda to solve for old builds.

## What's verified vs. not

Verified statically (no conda in the authoring environment): `meta.yaml` renders
+ parses; every dependency exists on the stated channel; the package imports and
the CLI runs on Python 3.14; no other removed-stdlib modules are used. **Not**
verified here: the actual `conda build` solve and the recipe `test:` block — run
those on a machine with conda.
