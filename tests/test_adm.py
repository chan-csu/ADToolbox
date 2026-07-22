from adtoolbox import __version__, adm


def test_version():
    assert __version__ == "1.1.0"


def test_adm_module_exports_model_class():
    assert hasattr(adm, "Model")
