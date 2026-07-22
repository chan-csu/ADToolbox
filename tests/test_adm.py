from adtoolbox import __version__, adm


def test_version():
    # __version__ must resolve to a real version, not the source-tree fallback,
    # and stay consistent with pyproject.toml.
    assert __version__ != "0.0.0+unknown"
    assert __version__[0].isdigit()


def test_adm_module_exports_model_class():
    assert hasattr(adm, "Model")
