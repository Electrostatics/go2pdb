# Retrieve PDB entries annotated with specific GO codes

This code retrieves [Protein DataBank](https://www.rcsb.org/) protein structure entries for proteins that have been annotated in [UniProt](https://www.uniprot.org/) with specific [Gene Ontology](http://geneontology.org/) codes.
This is useful for performing structural bioinformatics analysess (e.g., electrostatics comparisons, etc.) across proteins with similar functions.

This code also uses the [Gene Ontology Annotation database](https://www.ebi.ac.uk/GOA/) for enhanced annotation beyond what can be queried from UniProt.
This index is updated approximately every 4 weeks and can be [downloaded via FTP](ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/PDB/).
However, the database appears to have holes in it, so this code also searches UniProt for entries with matching GO codes that also have PDB structures.

## Installation and dependencies

After creating and activating a virtual environment (e.g., with [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://virtualenv.pypa.io/en/latest/)), you can install the code and most of its dependencies by running

```bash
pip install .
```

from the top of the source directory.
The code assumes that you have the docker program available in your path (i.e., can be run from the command line).

More information about the code use can be obtained by running

```bash
python -m go2pdb --help
```

## Example use

For more information about any of the commands below, using the `--help` option.
For example:

```bash
python -m go2pdb --help
python -m go2pdb search --help
```

### Finding structures with specific functions

For example, suppose we wanted to find the structures of all "nickel-binding proteins".
We will define this set of proteins as those gene products that have the Gene Ontology annotation [GO:0016151 "nickel cation binding"](https://www.ebi.ac.uk/QuickGO/term/GO:0016151).
However, if we're not completely sure that this GO code covers all relevant proteins, we can also add a PDB keyword search (e.g., for "nickel").

```bash
python -m go2pdb search --pdb-keyword NICKEL --search-goa GO:0016151
```

The `--search-goa` option adds search results from GOA.
Even though GOA appears to be incomplete, it provides useful additional information in its results (when present).

By default, the command above will produce a `search-output.xlsx` file that includes the results of your search.

### Comparing results by sequence similarity and identity

For most analyses, it is useful to start by grouping structures with similar sequences.
Sequence alignment is performed using the [BLAST Docker container](https://hub.docker.com/r/ncbi/blast).

**NOTE**: To perform this sequence analysis, you need to have [Docker](https://www.docker.com/) (e.g., [Docker Desktop](https://www.docker.com/products/docker-desktop)) installed on your computer.

The first step is to run BLAST on the existing sequences:

```bash
python -m go2pdb blast
```

This command consumes the `search-output.xlsx` file from the search step and produces a `blast-output.xlsx` file with pairwise matches between sequences.
Note that the results are filtered based on similarity and identity cutoffs; run with the `--help` option for more information.
The ``blast-output.xlsx` file can be used with graph visualization tools for qualitative insight into the relationships between PDB entries.
