# `optimize`

Parameter estimation and model validation. Four optimizer backends share one interface;
see the [parameter tuning guide](Optimization.md) for how to choose between them.

```python
from adtoolbox import optimize
```

!!! info "Optional dependencies"
    `ScipyOptimizer` works out of the box. The others need extras:
    `pip install "adtoolbox[blackbox]"` (OpenBox), `"adtoolbox[genetic]"` (PyGAD), or
    `"adtoolbox[surrogate]"` (PyTorch). Install all of them with `"adtoolbox[optimize]"`.

---

## Search space

### ParameterSpec

::: adtoolbox.optimize.ParameterSpec

### OptimizationRecord

::: adtoolbox.optimize.OptimizationRecord

---

## Base interface

Shared machinery: input validation, parameter-name resolution, model preparation,
objective evaluation, history bookkeeping, and JSON persistence.

::: adtoolbox.optimize.Optimizer

---

## Optimizer backends

### ScipyOptimizer

::: adtoolbox.optimize.ScipyOptimizer

### BlackBoxOptimizer

::: adtoolbox.optimize.BlackBoxOptimizer

### GeneticOptimizer

::: adtoolbox.optimize.GeneticOptimizer

### SurrogateOptimizer

::: adtoolbox.optimize.SurrogateOptimizer

---

## Validation

### validate_model

::: adtoolbox.optimize.validate_model

### calculate_fit_stats

::: adtoolbox.optimize.calculate_fit_stats
