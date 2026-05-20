# Reference Data

This folder contains clean reference JSON files that show the preferred ADToolbox database layout.

## Files

| File | Purpose | Top-level keys |
| --- | --- | --- |
| `models.json` | ADM model definitions. Each model entry contains `species`, `reactions`, `model_parameters`, `base_parameters`, `initial_conditions`, and `inlet_conditions`. | `adm1`, `e_adm` |
| `feeds.json` | Feed definitions keyed by feed name. Converted from `/Users/parsaghadermarzi/Desktop/Academics/Projects/Database/ADToolbox/feed_db.tsv`. Each entry matches the `core.Feed` constructor. | feed names |
| `experiments.json` | Experiment definitions keyed by experiment name. Converted from `/Users/parsaghadermarzi/Desktop/Academics/Projects/Database/ADToolbox/experimental_data_references.json`. Each entry matches the `core.Experiment` constructor shape. | experiment names |

The rule is one file per database type, with the first JSON layer used as the stable identifier for each item.

Example:

```json
{
  "adm1": {
    "species": [],
    "reactions": [],
    "model_parameters": {},
    "base_parameters": {},
    "initial_conditions": {},
    "inlet_conditions": {}
  }
}
```
