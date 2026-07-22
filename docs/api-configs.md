# `configs`

Configuration objects that tell every other module where files live and which thresholds
to use. Each config derives its paths from the directory you pass in, so the same code can
run against a scratch directory, a shared database, or a cluster filesystem without any
global state.

```python
from adtoolbox import configs
```

See the [API overview](API.md) for how configs pair with `core` objects.

## Helpers

::: adtoolbox.configs.adm_parameter_paths

## Database

Paths and remote URLs for reaction, compound, protein, feed, and study databases, plus the
ADM parameter bundle.

::: adtoolbox.configs.Database

## Metagenomics

Directories, alignment thresholds, and container images used by the metagenomics pipeline.
Defaults are derived from a `Database` config unless you override them.

::: adtoolbox.configs.Metagenomics

## Annotation

::: adtoolbox.configs.Annotation

## Documentation

::: adtoolbox.configs.Documentation

## Studies

::: adtoolbox.configs.Studies

## Utils

Container images and Slurm defaults used by the helpers in
[`utils`](api-utils.md).

::: adtoolbox.configs.Utils
