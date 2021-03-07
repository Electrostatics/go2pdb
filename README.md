# Retrieve PDB entries annotated with specific GO codes

This code retrieves [Protein DataBank](https://www.rcsb.org/) protein structure entries for proteins that have been annotated in the [Gene Ontology Annotation database](https://www.ebi.ac.uk/GOA/) with specific [Gene Ontology](http://geneontology.org/) codes.
This index is updated approximately every 4 weeks and can be [downloaded via FTP](ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/PDB/).

## Example use

For example, suppose we wanted to find the structures of all "toxic proteins".
We will define this set of proteins as those gene products that have the Gene Ontology annotation [GO:0090729 "toxin activity"](http://www.informatics.jax.org/vocab/gene_ontology/GO:0090729).

After installation (see below), the code can be run with
```
python -m pipeline
```

which executes the following sequence of steps:

1. Download an updated copy of ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/PDB/goa_pdb.gaf.gz.
2. This is a tab-delimited file with the following columns of interest:  `(2)` PDB ID, `(4)` qualifier on annotation, `(5)` GO annotation.  We want to extract column `(2)` for all entries where `(5)` is "GO:0090729" and `(4)` is not `NOT`.
3. Retrieve the protein structure from (for example) https://files.rcsb.org/download/1FAS.pdb and the protein sequence from (for example) https://www.rcsb.org/fasta/entry/1FAS.
4. Aggregate all of the FASTA sequences into a single file `all-toxins.fasta`.
If searches have already been performed for some sequences, you may want to create a second file `new-toxins.fasta` that contains only the sequences that haven't been searched.
5. Compile `toxins.fasta` into a BLAST database with
```
makeblastdb -in toxins.fasta -dbtype prot -out toxins_db
```
6. Search for sequence matches (`##` should be set to your target expectation value cutoff; probably 0):
```
blastp -db toxins_db -query new-toxins.fasta -max_hsps 1 -evalue ## -outfmt 5 -out output.xml
```
7. Use [BioPython](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec125) to parse the BLAST output.
8. Generate and store clusters of sequences based on sequence similarity.

The results are stored in `results.xlsx`.

More information about the code use can be obtained by running
```
python -m pipeline --help
```

## Installation and dependencies

After creating and activating a virtual environment (e.g., with conda or virtualenv), you can install the code and most of its dependencies by running
```
pip install .
```
from the top of the source directory.
The code assumes that you have the docker program available in your path (i.e., can be run from the command line).
