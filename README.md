# Retrieve PDB entries annotated with specific GO codes

This code retrieves [Protein DataBank](https://www.rcsb.org/) protein structure entries for proteins that have been annotated in [UniProt](https://www.uniprot.org/) with specific [Gene Ontology](http://geneontology.org/) codes.
This is useful for performing structural bioinformatics analyses (e.g., electrostatics comparisons, etc.) across proteins with similar functions.

This code also uses the [Gene Ontology Annotation database](https://www.ebi.ac.uk/GOA/) for enhanced annotation beyond what can be queried from UniProt.
This index is updated approximately every 4 weeks and can be [downloaded via FTP](ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/PDB/).
However, the database appears to have holes in it, so this code also searches UniProt for entries with matching GO codes that also have PDB structures.

## Support

This software is supported by the National Institutes of Health (grant GM069702).

## Installation and dependencies

After creating and activating a virtual environment (e.g., with [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://virtualenv.pypa.io/en/latest/)), you can install the code and most of its dependencies by running

```bash
pip install .
```

from the top of the source directory.
The code assumes that you have the docker program available in your path (i.e., can be run from the command line).

More information about the code use can be obtained by running

```bash
go2pdb --help
```

## Example use

For more information about any of the commands below, using the `--help` option.
For example:

```bash
go2pdb --help
go2pdb search --help
```

### Finding structures with specific functions

For example, suppose we wanted to find the structures of all "nickel-binding proteins".
We will define this set of proteins as those gene products that have the Gene Ontology annotation [GO:0016151 "nickel cation binding"](https://www.ebi.ac.uk/QuickGO/term/GO:0016151).
However, if we're not completely sure that this GO code covers all relevant proteins, we can also add a PDB keyword search (e.g., for "nickel").

```bash
go2pdb search --pdb-keyword NICKEL --search-goa GO:0016151
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
go2pdb blast
```

This command consumes the `search-output.xlsx` file from the search step and produces a `blast-output.xlsx` file with pairwise matches between sequences.
Note that the results are filtered based on similarity and identity cutoffs; run with the `--help` option for more information.
The ``blast-output.xlsx` file can be used with graph visualization tools for qualitative insight into the relationships between PDB entries.

### Clustering and summarizing results

Running

```bash
go2pdb cluster
```

will cluster sequences based on sequence identity.
The metric used for clustering can be changed with the `--cluster-metric` option and the cutoff for clustering can be changed with the `--metric-cutoff` option.
This clustering step will produce a simple table in `cluster-output.xlsx` that associates PDB chains with sequence-based clusters.
The sequence-based clusters are named by a representative protein in each cluster.

The cluster information can be merged with the search results by running

```bash
go2pdb summarize
```

which will produce a joined table in `summary-output.xlsx`.

The default behavior for both of these commands can be modified with options described in the `--help` option output.

The spreadsheet output contains the following columns:

Columns | Description
------- | -----------
PDB ID, PDB description, PDB title | Basic information about the structure from the Protein Data Bank (PDB)
PDB deposit date | The date the structure was added to the PDB (sort by this for the newest structures)
PDB method, PDB resolution (A) | Information about the experimental method for determining the structure; sort by resolution for the highest-refined structures
PDB chain ID, PDB strand ID(s), PDB strand type, PDB strand sequence | Information about a specific strand (subunit) of the protein.  This is the sequence used for BLAST comparisons between proteins.
PDB keyword match | If a PDB keyword search was performed, this gives the matching keyword (e.g., "NICKEL") used for the search.  When this is blank, it means that "NICKEL" did not occur in the PDB keywords for this structure.
UniProt entry ID, UniProt entry name, UniProt protein names | Description of the UniProt entry for this strand.  Often the UniProt description is more informative than the PDB description and the UniProt online entry contains links to many useful tools for analyzing the sequence.  When these fields are blank, it means the structure was not identified from a GO code match (i.e., was found from a PDB keyword search instead).
UniProt GO code | Matching GO code from search.
GOA qualifiers, GOA GO code, GOA DB reference, GOA evidence, GOA additional evidence, GOA taxon ID, GOA annotation date, GOA assigned by | If a search of the Gene Ontology Annotation (GOA) database was performed, this provides additional information about how the GO assignment was made.  The GOA database appears to be incomplete so blank entries should not be interpreted as lack of evidence for the GO annotation.
Cluster representative, Cluster description | These are the clusters assigned to each strand based on BLAST sequence identity with a 90% identity threshold.  Sorting by this column is useful for focusing on a specific type of protein.
