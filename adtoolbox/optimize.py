from __future__ import annotations

from abc import ABC, abstractmethod
from collections import namedtuple
from dataclasses import asdict, dataclass, field
import importlib
import json
import logging
import pathlib
from typing import Any, Iterable, Literal, Mapping, Sequence

from adtoolbox import adm, core
import numpy as np
import plotly
import plotly.express as px
import plotly.graph_objects as go
import polars as pl


Validation = namedtuple("Validation", ("r_squared", "rmse"))
ParameterTarget = Literal[
    "auto",
    "model_parameters",
    "base_parameters",
    "initial_conditions",
    "inlet_conditions",
]


def _require_package(module_name: str, extra: str):
    module = importlib.util.find_spec(module_name)
    if module is None:
        raise ImportError(
            f"This optimizer requires `{module_name}`. Install it with "
            f"`pip install adtoolbox[{extra}]`."
        )
    return importlib.import_module(module_name)


@dataclass(frozen=True)
class ParameterSpec:
    """Bounds and optional default value for one optimized parameter."""

    name: str
    lower: float
    upper: float
    default: float | None = None

    @classmethod
    def from_value(cls, name: str, value: "ParameterSpec | Sequence[float] | Mapping[str, float]") -> "ParameterSpec":
        if isinstance(value, ParameterSpec):
            return value
        if isinstance(value, Mapping):
            return cls(
                name=name,
                lower=float(value["lower"]),
                upper=float(value["upper"]),
                default=None if value.get("default") is None else float(value["default"]),
            )
        if len(value) not in (2, 3):
            raise ValueError(f"Search-space entry for `{name}` must be (lower, upper) or (lower, upper, default).")
        lower, upper = float(value[0]), float(value[1])
        default = None if len(value) == 2 else float(value[2])
        return cls(name=name, lower=lower, upper=upper, default=default)

    def __post_init__(self) -> None:
        if self.upper <= self.lower:
            raise ValueError(f"Upper bound for `{self.name}` must be greater than the lower bound.")
        if self.default is not None and not self.lower <= self.default <= self.upper:
            raise ValueError(f"Default value for `{self.name}` must be inside the search-space bounds.")


@dataclass
class OptimizationRecord:
    step: int
    parameters: dict[str, float]
    cost: float
    metadata: dict[str, Any] = field(default_factory=dict)


class Optimizer(ABC):
    """Common interface for ADToolbox parameter optimizers."""

    #: Cost assigned to a candidate that cannot be simulated (physically
    #: infeasible parameters or a non-finite objective). Large enough to be
    #: rejected by any population-based search without producing overflow.
    _INFEASIBLE_PENALTY = 1e30

    def __init__(
        self,
        base_model: adm.Model,
        train_data: Iterable[core.Experiment],
        search_space: Mapping[str, ParameterSpec | Sequence[float] | Mapping[str, float]],
        *,
        parameter_target: ParameterTarget = "auto",
        fitness_mode: str = "sum_squared_error",
        ode_method: str = "BDF",
        random_state: int | None = None,
        logger: logging.Logger | None = None,
    ) -> None:
        self.base_model = base_model
        self.train_data = list(train_data)
        self.search_space = {
            name: ParameterSpec.from_value(name, value)
            for name, value in search_space.items()
        }
        self.parameter_target = parameter_target
        self.fitness_mode = fitness_mode
        self.ode_method = ode_method
        self.random_state = random_state
        self.rng = np.random.default_rng(random_state)
        self.logger = logger or logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        self.history: list[OptimizationRecord] = []
        self.optimized_parameters: dict[str, float] | None = None
        self.optimized_model: adm.Model | None = None
        self.best_cost: float | None = None

        self._validate_inputs()

    @property
    def parameter_names(self) -> list[str]:
        """Optimized parameter names, in the order used by parameter vectors."""
        return list(self.search_space)

    @property
    def bounds(self) -> np.ndarray:
        """Search-space bounds as an ``(n_parameters, 2)`` array of lower/upper pairs."""
        return np.array(
            [[spec.lower, spec.upper] for spec in self.search_space.values()],
            dtype=float,
        )

    @property
    def best_record(self) -> OptimizationRecord | None:
        """The lowest-cost record seen so far, or None if nothing has been evaluated."""
        if not self.history:
            return None
        return min(self.history, key=lambda record: record.cost)

    def _validate_inputs(self) -> None:
        if not self.train_data:
            raise ValueError("At least one experiment is required for optimization.")
        if not self.search_space:
            raise ValueError("Search space cannot be empty.")
        if self.fitness_mode != "sum_squared_error":
            raise ValueError("Only `sum_squared_error` fitness is currently supported.")
        if self.parameter_target == "auto":
            missing = [
                name for name in self.search_space
                if self._parameter_location(self.base_model, name) is None
            ]
        else:
            available = self._target_keys(self.base_model, self.parameter_target)
            missing = [name for name in self.search_space if name not in available]
        if missing:
            joined = ", ".join(missing)
            raise ValueError(f"Search-space parameter(s) not found in model: {joined}")

    def _target_keys(self, model: adm.Model, target: ParameterTarget) -> set[str]:
        if target == "model_parameters":
            return set(model.model_parameters)
        if target == "base_parameters":
            return set(model.base_parameters)
        if target == "initial_conditions":
            return set(getattr(model, "_ic", {}))
        if target == "inlet_conditions":
            return set(getattr(model, "_inc", {})) | set(model.species)
        return set()

    def _parameter_location(self, model: adm.Model, name: str) -> ParameterTarget | None:
        if name in model.model_parameters:
            return "model_parameters"
        if name in model.base_parameters:
            return "base_parameters"
        if name in getattr(model, "_ic", {}):
            return "initial_conditions"
        if name in getattr(model, "_inc", {}) or name in model.species:
            return "inlet_conditions"
        return None

    def _coerce_parameters(self, parameters: Mapping[str, float] | Sequence[float] | np.ndarray) -> dict[str, float]:
        if isinstance(parameters, Mapping):
            missing = set(self.search_space) - set(parameters)
            unknown = set(parameters) - set(self.search_space)
            if missing or unknown:
                problems = []
                if missing:
                    problems.append(f"missing: {', '.join(sorted(missing))}")
                if unknown:
                    problems.append(f"unknown: {', '.join(sorted(unknown))}")
                raise ValueError("Invalid parameter set (" + "; ".join(problems) + ").")
            return {name: float(parameters[name]) for name in self.parameter_names}

        vector = np.asarray(parameters, dtype=float).reshape(-1)
        if vector.size != len(self.search_space):
            raise ValueError(f"Expected {len(self.search_space)} parameters, got {vector.size}.")
        return dict(zip(self.parameter_names, vector.tolist()))

    def parameters_to_vector(self, parameters: Mapping[str, float]) -> np.ndarray:
        """Convert a parameter mapping into a vector ordered by `parameter_names`."""
        parameters = self._coerce_parameters(parameters)
        return np.array([parameters[name] for name in self.parameter_names], dtype=float)

    def vector_to_parameters(self, vector: Sequence[float] | np.ndarray) -> dict[str, float]:
        """Convert a parameter vector back into a name-to-value mapping."""
        return self._coerce_parameters(vector)

    def clip_vector(self, vector: Sequence[float] | np.ndarray) -> np.ndarray:
        """Clip a parameter vector element-wise into the search-space bounds."""
        values = np.asarray(vector, dtype=float).reshape(-1)
        return np.clip(values, self.bounds[:, 0], self.bounds[:, 1])

    def default_vector(self) -> np.ndarray:
        """Parameter vector of each spec's default, or its bound midpoint if unset."""
        return np.array(
            [
                spec.default if spec.default is not None else (spec.lower + spec.upper) / 2
                for spec in self.search_space.values()
            ],
            dtype=float,
        )

    def random_vector(self) -> np.ndarray:
        """Draw a uniformly random parameter vector from within the bounds."""
        return self.rng.uniform(self.bounds[:, 0], self.bounds[:, 1])

    def prepare_model(
        self,
        parameters: Mapping[str, float] | Sequence[float] | np.ndarray,
        experiment: core.Experiment | None = None,
    ) -> adm.Model:
        """Copy the base model and apply candidate parameters.

        When an experiment is given, its feed, base parameters, and initial
        concentrations are applied as well, the initial value of every measured
        variable is set to its observed value at time zero, and any state listed
        in the experiment's `constants` is pinned as a control state.

        Args:
            parameters: Candidate values, as a mapping or a vector.
            experiment: Optional experiment whose conditions should be applied.

        Returns:
            adm.Model: A new model instance ready to solve. The base model is
                never mutated.
        """
        model = self.base_model.copy()
        self._apply_parameters(model, self._coerce_parameters(parameters))
        if experiment is not None:
            if experiment.base_parameters:
                model.update_parameters(base_parameters=experiment.base_parameters)
            model.feed = experiment.feed
            ic = self._experiment_initial_conditions(experiment)
            model.update_parameters(initial_conditions=ic)
            model.control_state = {
                key: ic[key]
                for key in experiment.constants
                if key in ic
            }
        return model

    def _experiment_initial_conditions(self, experiment: core.Experiment) -> dict[str, float]:
        initial = dict(getattr(self.base_model, "_ic", {}))
        initial.update(experiment.initial_concentrations)
        for index, variable in enumerate(experiment.variables):
            initial[variable] = float(experiment.data[0, index])
        return initial

    def _apply_parameters(self, model: adm.Model, parameters: Mapping[str, float]) -> None:
        updates: dict[str, dict[str, float]] = {
            "model_parameters": {},
            "base_parameters": {},
            "initial_conditions": {},
            "inlet_conditions": {},
        }
        for name, value in parameters.items():
            target = self.parameter_target
            if target == "auto":
                target = self._parameter_location(model, name)
            if target is None:
                raise ValueError(f"Could not determine target for parameter `{name}`.")
            if target == "inlet_conditions" and name.endswith("_in"):
                name = name[:-3]
            updates[target][name] = float(value)

        clean_updates = {key: value for key, value in updates.items() if value}
        if clean_updates:
            model.update_parameters(**clean_updates)

    def evaluate(self, parameters: Mapping[str, float] | Sequence[float] | np.ndarray) -> float:
        """Score one parameter set against every training experiment.

        Each experiment is simulated at its own measurement time points and
        compared to the observed data. Calling this directly does not add to
        the optimizer history.

        Args:
            parameters: Candidate values, as a mapping or a vector.

        Returns:
            float: Sum of squared residuals across all experiments, time
                points, and measured variables. Lower is better.
        """
        parameters = self._coerce_parameters(parameters)
        total = 0.0
        for experiment in self.train_data:
            try:
                model = self.prepare_model(parameters, experiment)
                rows = [model.species.index(variable) for variable in experiment.variables]
                solution = model.solve_model(np.array(experiment.time), method=self.ode_method)
                prediction = np.asarray(solution.y[rows, :], dtype=float).T
                residual = prediction - experiment.data
                cost = float(np.sum(np.square(residual)))
            except Exception as exc:
                # A candidate the optimizer proposes can be physically
                # infeasible (e.g. fermentation fractions that make a derived
                # coefficient negative). That is not a bug in the model — it is
                # a point outside the feasible region — so we penalise it
                # heavily instead of letting it abort the whole search.
                self.logger.debug(
                    "penalising infeasible candidate on %s: %s", experiment.name, exc
                )
                cost = self._INFEASIBLE_PENALTY
            if not np.isfinite(cost):
                cost = self._INFEASIBLE_PENALTY
            total += cost
        return total

    def _evaluate_and_record(
        self,
        parameters: Mapping[str, float] | Sequence[float] | np.ndarray,
        *,
        metadata: Mapping[str, Any] | None = None,
    ) -> float:
        parameters = self._coerce_parameters(parameters)
        cost = self.evaluate(parameters)
        self.record(parameters, cost, metadata=metadata)
        return cost

    def record(
        self,
        parameters: Mapping[str, float] | Sequence[float] | np.ndarray,
        cost: float,
        *,
        metadata: Mapping[str, Any] | None = None,
    ) -> OptimizationRecord:
        """Append an evaluated point to the history.

        If the cost improves on the current best, `best_cost`,
        `optimized_parameters`, and `optimized_model` are updated.

        Args:
            parameters: The evaluated values, as a mapping or a vector.
            cost: The objective value for those parameters.
            metadata: Optional free-form annotations, such as which backend
                produced the point.

        Returns:
            OptimizationRecord: The record that was appended.
        """
        record = OptimizationRecord(
            step=len(self.history),
            parameters=self._coerce_parameters(parameters),
            cost=float(cost),
            metadata=dict(metadata or {}),
        )
        self.history.append(record)
        if self.best_cost is None or record.cost < self.best_cost:
            self.best_cost = record.cost
            self.optimized_parameters = record.parameters.copy()
            self.optimized_model = self.prepare_model(record.parameters)
        return record

    def clear_history(self) -> None:
        """Discard all recorded evaluations and the current best result."""
        self.history.clear()
        self.optimized_parameters = None
        self.optimized_model = None
        self.best_cost = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the search space, best result, and full history to a dict."""
        return {
            "optimizer": self.__class__.__name__,
            "parameter_target": self.parameter_target,
            "fitness_mode": self.fitness_mode,
            "ode_method": self.ode_method,
            "search_space": {
                name: asdict(spec)
                for name, spec in self.search_space.items()
            },
            "optimized_parameters": self.optimized_parameters,
            "best_cost": self.best_cost,
            "history": [asdict(record) for record in self.history],
        }

    def save(self, path: str | pathlib.Path) -> pathlib.Path:
        """Write the optimizer state to a JSON file, creating parent directories.

        Args:
            path: Destination file path.

        Returns:
            pathlib.Path: The path that was written.
        """
        path = pathlib.Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2)
        return path

    def load(self, path: str | pathlib.Path) -> "Optimizer":
        """Restore optimizer state by replaying a saved history.

        The current history is cleared first, so `best_cost`,
        `optimized_parameters`, and `optimized_model` end up reflecting the
        saved run.

        Args:
            path: A JSON file previously written by `save`.

        Returns:
            Optimizer: This optimizer, to allow chaining.

        Raises:
            ValueError: If the saved parameter names do not match this
                optimizer's search space.
        """
        path = pathlib.Path(path)
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        saved_names = list(payload.get("search_space", {}))
        if saved_names and saved_names != self.parameter_names:
            raise ValueError("Saved optimizer state does not match this optimizer search space.")
        self.clear_history()
        for item in payload.get("history", []):
            self.record(item["parameters"], item["cost"], metadata=item.get("metadata"))
        return self

    @abstractmethod
    def optimize(self, **kwargs) -> Any:
        """Run the optimizer and update history, optimized_parameters, and optimized_model."""


class ScipyOptimizer(Optimizer):
    """Differential-evolution optimizer using SciPy only."""

    def optimize(
        self,
        *,
        maxiter: int = 100,
        popsize: int = 15,
        polish: bool = True,
        workers: int = 1,
        **kwargs,
    ):
        from scipy.optimize import differential_evolution

        def objective(vector: np.ndarray) -> float:
            return self._evaluate_and_record(
                vector,
                metadata={"optimizer": "scipy_differential_evolution"},
            )

        result = differential_evolution(
            objective,
            bounds=[tuple(bound) for bound in self.bounds],
            maxiter=maxiter,
            popsize=popsize,
            polish=polish,
            seed=self.random_state,
            workers=workers,
            **kwargs,
        )
        self.record(
            result.x,
            float(result.fun),
            metadata={"optimizer": "scipy_differential_evolution", "final": True},
        )
        return result


class BlackBoxOptimizer(Optimizer):
    """OpenBox-backed optimizer with the shared ADToolbox optimizer interface."""

    def optimize(self, *, parallel: bool = False, **kwargs):
        openbox = _require_package("openbox", "blackbox")
        space_module = importlib.import_module("openbox.space")
        space = space_module.Space()
        space.add_variables(
            [
                space_module.Real(
                    spec.name,
                    spec.lower,
                    spec.upper,
                    default_value=spec.default if spec.default is not None else (spec.lower + spec.upper) / 2,
                )
                for spec in self.search_space.values()
            ]
        )

        def objective(config):
            parameters = config.get_dictionary() if hasattr(config, "get_dictionary") else dict(config)
            return self._evaluate_and_record(
                parameters,
                metadata={"optimizer": "openbox"},
            )

        optimizer_cls = openbox.ParallelOptimizer if parallel else openbox.Optimizer
        optimizer = optimizer_cls(objective, space, **kwargs)
        result = optimizer.run()
        record = self.best_record
        if record is not None:
            self.optimized_parameters = record.parameters.copy()
            self.best_cost = record.cost
            self.optimized_model = self.prepare_model(record.parameters)
        return result


class GeneticOptimizer(Optimizer):
    """PyGAD-backed genetic optimizer with the shared ADToolbox optimizer interface."""

    def optimize(self, *, num_generations: int = 50, sol_per_pop: int = 20, **kwargs):
        pygad = _require_package("pygad", "genetic")

        def fitness(ga_instance, solution, solution_idx) -> float:
            cost = self._evaluate_and_record(
                solution,
                metadata={"optimizer": "pygad", "solution_idx": int(solution_idx)},
            )
            return 1.0 / (cost + 1e-12)

        ga = pygad.GA(
            fitness_func=fitness,
            num_genes=len(self.search_space),
            gene_space=[
                {"low": spec.lower, "high": spec.upper}
                for spec in self.search_space.values()
            ],
            num_generations=num_generations,
            sol_per_pop=sol_per_pop,
            **kwargs,
        )
        ga.run()
        solution, fitness_value, solution_idx = ga.best_solution()
        cost = 1.0 / max(float(fitness_value), 1e-12)
        self.record(
            solution,
            cost,
            metadata={"optimizer": "pygad", "final": True, "solution_idx": int(solution_idx)},
        )
        return ga


class SurrogateOptimizer(Optimizer):
    """Small neural surrogate optimizer for expensive ODE model evaluations."""

    def __init__(
        self,
        *args,
        hidden_size: int = 30,
        hidden_layers: int = 4,
        learning_rate: float = 1e-3,
        input_learning_rate: float = 1e-3,
        train_epochs: int = 500,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        torch = _require_package("torch", "surrogate")
        layers: list[Any] = []
        width = len(self.search_space)
        for _ in range(hidden_layers):
            layers.append(torch.nn.Linear(width, hidden_size))
            layers.append(torch.nn.Tanh())
            width = hidden_size
        layers.append(torch.nn.Linear(width, 1))
        self._torch = torch
        self.network = torch.nn.Sequential(*layers)
        self.learning_rate = learning_rate
        self.input_learning_rate = input_learning_rate
        self.train_epochs = train_epochs

    def _history_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.history:
            raise ValueError("Surrogate optimizer needs at least one evaluated point.")
        x = np.array(
            [
                [record.parameters[name] for name in self.parameter_names]
                for record in self.history
            ],
            dtype=np.float32,
        )
        y = np.array([record.cost for record in self.history], dtype=np.float32).reshape(-1, 1)
        return x, y

    def _train_surrogate(self) -> float:
        torch = self._torch
        x, y = self._history_arrays()
        inputs = torch.tensor(x, dtype=torch.float32)
        labels = torch.tensor(y, dtype=torch.float32)
        labels = torch.clamp(labels, max=1e6)
        optimizer = torch.optim.Adam(self.network.parameters(), lr=self.learning_rate)
        loss = torch.tensor(float("nan"))
        for _ in range(self.train_epochs):
            predicted = self.network(inputs)
            loss = torch.mean(torch.square(predicted - labels))
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
        return float(loss.detach().numpy())

    def _suggest_vector(self, *, grad_steps: int) -> np.ndarray:
        torch = self._torch
        start = self.parameters_to_vector(self.optimized_parameters) if self.optimized_parameters else self.default_vector()
        vector = torch.tensor(start, dtype=torch.float32, requires_grad=True)
        optimizer = torch.optim.Adam([vector], lr=self.input_learning_rate)
        lower = torch.tensor(self.bounds[:, 0], dtype=torch.float32)
        upper = torch.tensor(self.bounds[:, 1], dtype=torch.float32)
        for _ in range(grad_steps):
            loss = torch.mean(self.network(vector))
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            with torch.no_grad():
                vector.copy_(torch.maximum(torch.minimum(vector, upper), lower))
        return vector.detach().numpy()

    def optimize(
        self,
        *,
        n_steps: int = 100,
        initial_points: int = 8,
        grad_steps: int = 20,
        perturbation_scale: float = 0.05,
        save_every: int | None = None,
        history_path: str | pathlib.Path | None = None,
    ) -> list[OptimizationRecord]:
        if not self.history:
            for index in range(initial_points):
                self._evaluate_and_record(
                    self.random_vector(),
                    metadata={"optimizer": "surrogate", "phase": "initial", "index": index},
                )

        for step in range(n_steps):
            loss = self._train_surrogate()
            candidate = self.clip_vector(self._suggest_vector(grad_steps=grad_steps))
            if self.history and len(self.history) > 1:
                recent = [record.cost for record in self.history[-2:]]
                if abs(recent[-1] - recent[-2]) < 1e-9:
                    spread = np.maximum(np.abs(candidate) * perturbation_scale, perturbation_scale)
                    candidate = self.clip_vector(self.rng.normal(candidate, spread))
            cost = self._evaluate_and_record(
                candidate,
                metadata={"optimizer": "surrogate", "phase": "search", "step": step, "training_loss": loss},
            )
            self.logger.info("Surrogate step %s/%s cost=%s best=%s", step + 1, n_steps, cost, self.best_cost)
            if history_path is not None and save_every and (step + 1) % save_every == 0:
                self.save(history_path)

        if history_path is not None:
            self.save(history_path)
        return self.history


def validate_model(
    model: adm.Model,
    data: core.Experiment | Iterable[core.Experiment],
    plot: bool = False,
    show_extra_states: Iterable[str] | None = None,
    ode_solver: str = "Radau",
) -> tuple[dict[str, pl.DataFrame], plotly.graph_objs.Figure | None]:
    """
    Compare model predictions against one or more experiments.
    """
    palette = px.colors.qualitative.Plotly
    fig = None

    if not isinstance(data, core.Experiment):
        experiments = list(data)
        if not experiments:
            raise ValueError("At least one experiment is required.")
        fig = go.Figure()
        first = experiments[0]
        ic = pl.DataFrame([experiment.initial_concentrations for experiment in experiments]).mean().to_dicts()[0]
        ic.update({variable: first.data[0, idx] for idx, variable in enumerate(first.variables)})
        model.control_state = {key: ic[key] for key in first.constants if key in ic}
        model.update_parameters(base_parameters=first.base_parameters)
        model.update_parameters(initial_conditions=ic)
        all_time_points = np.array(sorted(set(sum([experiment.time for experiment in experiments], start=[]))))
        solution = model.solve_model(all_time_points, method=ode_solver)
        out = {
            "model": pl.DataFrame(
                solution.y[[model.species.index(variable) for variable in first.variables], :].T,
                schema=first.variables,
            ).with_columns(pl.Series("time", all_time_points.tolist()))
        }
        for idx, variable in enumerate(first.variables):
            fig.add_trace(
                go.Scatter(
                    x=out["model"]["time"],
                    y=out["model"][variable],
                    name=variable,
                    mode="lines",
                    line=dict(color=palette[idx]),
                )
            )
        rows = []
        for experiment in experiments:
            for row_index, time_point in enumerate(experiment.time):
                for col_index, variable in enumerate(experiment.variables):
                    rows.append(
                        {"time": time_point, "variable": variable, "value": experiment.data[row_index, col_index]}
                    )
        observed = pl.DataFrame(rows)
        comps = observed["variable"].unique().to_list()
        times = observed["time"].unique().to_list()
        for idx, comp in enumerate(comps):
            t, y, e = [], [], []
            for time_point in times:
                group = observed.filter((pl.col("time") == time_point) & (pl.col("variable") == comp))
                t.append(time_point)
                y.append(group["value"].mean())
                e.append((group["value"].std() or 0.0) / np.sqrt(group.height))
            fig.add_trace(
                go.Scatter(
                    x=t,
                    y=y,
                    error_y=dict(type="data", array=e, color=palette[idx], visible=True),
                    mode="markers",
                    marker=dict(color=palette[idx]),
                    name=f"{comp} Observed",
                )
            )
        if show_extra_states:
            for idx, extra in enumerate(show_extra_states):
                fig.add_trace(
                    go.Scatter(
                        x=out["model"]["time"],
                        y=solution.y[model.species.index(extra)],
                        name=extra,
                        mode="lines",
                        line=dict(color=palette[idx + len(first.variables)]),
                    )
                )
        if plot:
            fig.show(renderer="svg")
        return out, fig

    ic = data.initial_concentrations.copy()
    ic.update({variable: data.data[0, idx] for idx, variable in enumerate(data.variables)})
    model.control_state = {key: ic[key] for key in data.constants if key in ic}
    model.update_parameters(base_parameters=data.base_parameters)
    model.update_parameters(initial_conditions=ic)
    solution = model.solve_model(np.array(data.time), method=ode_solver)
    out = {
        "model": pl.DataFrame(
            solution.y[[model.species.index(variable) for variable in data.variables], :].T,
            schema=data.variables,
        ).with_columns(pl.Series("time", np.array(data.time).tolist())),
        "data": pl.DataFrame(data.data, schema=data.variables).with_columns(pl.Series("time", data.time)),
    }
    if plot:
        fig = go.Figure()
        for idx, variable in enumerate(data.variables):
            fig.add_trace(
                go.Scatter(
                    x=out["model"]["time"],
                    y=out["model"][variable],
                    name=variable,
                    mode="lines",
                    line=dict(color=palette[idx]),
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=out["data"]["time"],
                    y=out["data"][variable],
                    name=f"{variable} observed",
                    mode="markers",
                    marker=dict(color=palette[idx]),
                )
            )
        if show_extra_states:
            for idx, extra in enumerate(show_extra_states):
                fig.add_trace(
                    go.Scatter(
                        x=out["model"]["time"],
                        y=solution.y[model.species.index(extra)],
                        name=extra,
                        mode="lines",
                        line=dict(color=palette[idx + len(data.variables)]),
                    )
                )
        fig.show(renderer="svg")

    return out, fig


def calculate_fit_stats(model: adm.Model, data: Iterable[core.Experiment]) -> Validation:
    """Calculate RMSE and R-squared metrics for a set of experiments."""
    predicted_values = []
    observed_values = []
    for experiment in data:
        formatted_data = validate_model(model, experiment)[0]
        model_frame, data_frame = formatted_data["model"], formatted_data["data"]
        for column in model_frame.columns:
            if column == "time":
                continue
            predicted_values.extend(model_frame[column])
            observed_values.extend(data_frame[column])
    predicted = np.array(predicted_values)
    observed = np.array(observed_values)
    residual_sum = np.sum(np.square(observed - predicted))
    total_sum = np.sum(np.square(observed - np.mean(observed)))
    return Validation(
        r_squared=1 - (residual_sum / total_sum),
        rmse=np.sqrt(residual_sum),
    )
