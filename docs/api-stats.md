# `stats`

Small statistics helpers for feature tables produced by the metagenomics pipeline.

```python
from adtoolbox import stats
```

Both functions operate on Polars DataFrames and take a `feature_axis` argument: use `1`
when features are columns, `0` when features are rows.

---

## calculate_dist

Pairwise distance matrix, computed in parallel across a thread pool.

::: adtoolbox.stats.calculate_dist

## scaler

::: adtoolbox.stats.scaler

## Pair

::: adtoolbox.stats.Pair
