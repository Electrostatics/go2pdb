"""Process PDB entries with specific Gene Ontology annotations."""
import pkg_resources
from importlib import metadata
from pathlib import Path
from configparser import ConfigParser, ExtendedInterpolation


__version__ = metadata.version("go_pdb")


# Default command line options
DATA_DIR = Path("data")
CIF_DIR = DATA_DIR / Path("cif")
FASTA_DIR = DATA_DIR / Path("fasta")
GOA_DIR = DATA_DIR
PDB_DIR = DATA_DIR / Path("pdb")
BLAST_DIR = DATA_DIR / Path("blast")
