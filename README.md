# Retrieve PDB entries annotated with specific GO codes

This code retrieves [Protein DataBank](https://www.rcsb.org/) protein structure entries for proteins that have been annotated in the [Gene Ontology Annotation database](https://www.ebi.ac.uk/GOA/) with specific [Gene Ontology](http://geneontology.org/) codes.
This index is updated approximately every 4 weeks and can be [downloaded via FTP](ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/PDB/).
However, the database appears to have holes in it, so this code also searches UniProt for entries with matching GO codes that also have PDB structures.

## Example use

For example, suppose we wanted to find the structures of all "nickel-binding proteins".
We will define this set of proteins as those gene products that have the Gene Ontology annotation [GO:0016151 "nickel cation binding"](https://www.ebi.ac.uk/QuickGO/term/GO:0016151).
However, if we're not completely sure that this GO code covers all relevant proteins, we can also add a PDB keyword search (e.g., for "nickel").

After installation (see below), the search for matching structures can be run with
```
python -m go2pdb search --pdb-keyword NICKEL --search-goa GO:0016151
```
The `--search-goa` option adds search results from GOA.
Even though GOA appears to be incomplete, it provides useful additional information in its results (when present).

More information about the code use can be obtained by running
```
python -m go2pdb --help
```

## Installation and dependencies

After creating and activating a virtual environment (e.g., with `conda create` or `virtualenv`), you can install the code and most of its dependencies by running
```
pip install .
```
from the top of the source directory.
The code assumes that you have the docker program available in your path (i.e., can be run from the command line).
