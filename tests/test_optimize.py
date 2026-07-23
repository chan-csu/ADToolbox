import inspect

import numpy as np

from adtoolbox import __version__, core, optimize


class FakeSolution:
    def __init__(self, y):
        self.y = y


class FakeModel:
    species = ["S_x"]

    def __init__(self, k=0.0):
        self.model_parameters = {"k": k}
        self.base_parameters = {}
        self._ic = {"S_x": 0.0}
        self._inc = {}
        self.feed = None
        self.control_state = {}

    def copy(self):
        copied = FakeModel(self.model_parameters["k"])
        copied.base_parameters = self.base_parameters.copy()
        copied._ic = self._ic.copy()
        copied._inc = self._inc.copy()
        return copied

    def update_parameters(
        self,
        model_parameters=None,
        base_parameters=None,
        initial_conditions=None,
        inlet_conditions=None,
    ):
        if model_parameters:
            self.model_parameters.update(model_parameters)
        if base_parameters:
            self.base_parameters.update(base_parameters)
        if initial_conditions:
            self._ic.update(initial_conditions)
        if inlet_conditions:
            self._inc.update(inlet_conditions)

    def solve_model(self, time, method="BDF"):
        time = np.asarray(time, dtype=float)
        y = self._ic["S_x"] + self.model_parameters["k"] * time
        return FakeSolution(y.reshape(1, -1))


def _experiment():
    return core.Experiment(
        name="linear",
        time=[0, 1, 2],
        variables=["S_x"],
        data=[[1, 3, 5]],
        feed=core.Feed("test", carbohydrates=25, lipids=25, proteins=25, tss=50, si=25, xi=25),
    )


def test_version():
    assert __version__ != "0.0.0+unknown"
    assert __version__[0].isdigit()


def test_optimizer_is_abstract():
    assert inspect.isabstract(optimize.Optimizer)


def test_optimizer_evaluates_by_species_name():
    optimizer = optimize.ScipyOptimizer(
        base_model=FakeModel(),
        train_data=[_experiment()],
        search_space={"k": (0, 4)},
    )

    assert optimizer.evaluate({"k": 2.0}) == 0.0
    assert optimizer.evaluate([1.0]) > 0.0


def test_optimizer_tracks_and_loads_history(tmp_path):
    optimizer = optimize.ScipyOptimizer(
        base_model=FakeModel(),
        train_data=[_experiment()],
        search_space={"k": (0, 4)},
    )

    optimizer.record({"k": 1.0}, optimizer.evaluate({"k": 1.0}))
    optimizer.record({"k": 2.0}, optimizer.evaluate({"k": 2.0}))
    path = optimizer.save(tmp_path / "history.json")

    restored = optimize.ScipyOptimizer(
        base_model=FakeModel(),
        train_data=[_experiment()],
        search_space={"k": (0, 4)},
    ).load(path)

    assert restored.optimized_parameters == {"k": 2.0}
    assert restored.best_cost == 0.0
    assert len(restored.history) == 2
