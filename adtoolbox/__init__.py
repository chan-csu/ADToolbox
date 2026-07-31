import os
import sys

"""Project Setup for ADToolBox."""

from importlib.metadata import PackageNotFoundError, version as _pkg_version


def _resolve_version() -> str:
    # Single source of truth: pyproject.toml. Installed packages read it from
    # the recorded metadata; a source checkout reads pyproject.toml directly.
    try:
        return _pkg_version("adtoolbox")
    except PackageNotFoundError:
        pass
    try:
        import tomllib

        pyproject = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyproject.toml")
        with open(pyproject, "rb") as f:
            return tomllib.load(f)["tool"]["poetry"]["version"]
    except Exception:
        return "0.0.0+unknown"


__version__ = _resolve_version()
__all__=["adm","configs","__main__","cli","core","markers","optimize","pipeline","utils","PKG_DATA"]

sys.path.append(os.path.join(os.path.dirname(__file__)))

PKG_DATA=os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data")
