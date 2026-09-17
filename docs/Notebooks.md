# Example Notebooks

The `Examples/` directory in the repository holds executable notebooks for each major
workflow. Three of them are also rendered here as static walkthroughs.

## Available notebooks

<div class="grid cards" markdown>

-   :material-chart-bell-curve: **`modeling.ipynb`**

    ---

    Build a new ADM model or load an existing one, then solve, inspect, and visualize it.
    The best starting point if you care about simulation rather than sequencing.

    [:octicons-arrow-right-24: Rendered walkthrough](modeling.html)
    · [:octicons-mark-github-16: Source](https://github.com/chan-csu/ADToolbox/blob/main/Examples/modeling.ipynb)

-   :material-tune-variant: **`parameter_tuning.ipynb`**

    ---

    Reproduces the figures from the ADToolbox paper, including predicting AD products from
    16S data at time zero, and walks through the full parameter-tuning workflow.

    [:octicons-arrow-right-24: Rendered walkthrough](parameter_tuning.html)
    · [:octicons-mark-github-16: Source](https://github.com/chan-csu/ADToolbox/blob/main/Examples/parameter_tuning.ipynb)

-   :material-dna: **`metagenomics.ipynb`**

    ---

    End-to-end tutorial: from an SRA accession table to the NMDS + feature-enrichment
    figure. Runs the metagenomics pipeline (SRA → per-sample COD / marker profiles) via
    the Python API, then ordinates and compares the feedstock groups. Ends with the
    one-command CLI equivalent.

    [:octicons-arrow-right-24: Rendered walkthrough](metagenomics.html)
    · [:octicons-mark-github-16: Source](https://github.com/chan-csu/ADToolbox/blob/main/Examples/metagenomics.ipynb)

</div>

## Running them

=== "Locally"

    ```bash
    git clone https://github.com/chan-csu/ADToolbox.git
    cd ADToolbox
    pip install -e ".[optimize,dashboard]"
    pip install jupyterlab
    jupyter lab Examples/
    ```

=== "Binder"

    [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/chan-csu/ADToolbox/HEAD)

    No installation required. Escher map visualization is not available on Binder.

## Supporting data

The notebooks read from files that ship with the repository:

| Path | Contents |
| --- | --- |
| `Examples/ADM_Parameters/` | Individual e-ADM parameter JSON files. |
| `Examples/Studies/` | Metagenomics study metadata and experimental data references. |
| `Examples/feed_db.tsv` | A small feed database. |
| `Examples/toy_model.json`, `Examples/Toymap.json` | A minimal model and its Escher map. |
| `reference_data/models.json` | ADM1 and e-ADM parameter sets keyed by model name. |

!!! note "Notebooks that need databases"
    The metagenomics notebook expects the reference databases. Run
    `adtoolbox database download-all-databases --output-dir ./database` first, and point
    the notebook's `database_dir` at that directory. See the
    [Quickstart](Quickstart.md#2-download-the-reference-databases).
