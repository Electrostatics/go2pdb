"""Process PDB entries with specific Gene Ontology annotations."""
import pkg_resources
from importlib import metadata
from pathlib import Path
from configparser import ConfigParser, ExtendedInterpolation


__version__ = metadata.version("go2pdb")
