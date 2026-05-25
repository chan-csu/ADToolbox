import os
import sys

"""Project Setup for ADToolBox."""

__version__ = "1.0.0"
__all__=["adm","configs","__main__","cli","core","optimize","pipeline","utils","PKG_DATA"]

sys.path.append(os.path.join(os.path.dirname(__file__)))

PKG_DATA=os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data")
