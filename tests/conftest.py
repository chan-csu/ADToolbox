import importlib.util
import sys
import types


class _OptionalModule(types.ModuleType):
    def __init__(self, name):
        super().__init__(name)
        self.__file__ = f"<optional-stub-{name}>"
        self.__path__ = []

    def __getattr__(self, name):
        value = _OptionalModule(name)
        setattr(self, name, value)
        return value

    def __call__(self, *args, **kwargs):
        return _OptionalModule(self.__name__)


def _stub(module_name):
    sys.modules[module_name] = _OptionalModule(module_name)


if importlib.util.find_spec("dash") is None:
    for _module_name in ["dash", "dash.html", "dash.dcc", "dash.dash_table"]:
        _stub(_module_name)

if importlib.util.find_spec("dash_bootstrap_components") is None:
    _stub("dash_bootstrap_components")

if importlib.util.find_spec("dash_escher") is None:
    _stub("dash_escher")
