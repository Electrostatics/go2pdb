"""Manage BLAST calculations on sequences."""
import logging
import subprocess
import shutil
import json
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np


_LOGGER = logging.getLogger(__name__)
BLAST_TIMEOUT = 720000
BLAST_BUILD_CMD_FMT = (
    "docker run --rm -v {work_dir}:/work ncbi/blast "
    "makeblastdb -in /work/{sequence_file} -dbtype prot -parse_seqids "
    "-out /work/{db_file} -blastdb_version 5 -title 'sequences'"
)
BLAST_QUERY_CMD_FMT = (
    "docker run --rm -v {work_dir}:/work ncbi/blast "
    "blastp -query /work/{sequence_file} -max_hsps 1 -db /work/{db_file} "
    "-outfmt 5 -out /work/{output_file}"
)
OUTPUT_FILE = "output.xml"
SEQUENCE_FILE = "sequences.fasta"
DB_FILE = "blastdb"


def build_fasta(search_df) -> str:
    """Build a FASTA file from the sequences in search_df.

    :param pd.DataFrame search_df:  dataframe with sequences
    :returns:  string with FASTA data
    """
    fasta = ""
    search_df = search_df[["PDB chain ID", "PDB strand sequence",]]
    search_df = search_df.drop_duplicates()
    for _, row in search_df.iterrows():
        chain_id = row["PDB chain ID"]
        sequence = row["PDB strand sequence"]
        description = f"{chain_id}"
        _LOGGER.debug(f"Processing {description}...")
        try:
            seq = Seq(sequence)
        except TypeError:
            _LOGGER.warning(
                f"Failed to parse sequence {sequence} for {description} from "
                f"{chain_id}. This often happens when entries are withdrawn "
                f"from the PDB."
            )
            seq = None
        if seq is not None:
            fasta += SeqRecord(seq, id=description).format("fasta")
    return fasta


def build_db(sequence_file, blast_dir):
    """Build the BLAST database.

    :param str sequence_file:  path to FASTA file with sequences for DB
    :param Path blast_dir:  directory for BLAST results
    """
    # The temporary directory fixes problems with spaces in paths
    with TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        _LOGGER.debug(
            f"Using temporary directory {temp_dir} to build BLAST database."
        )
        shutil.copy(sequence_file.resolve(), Path(temp_dir) / SEQUENCE_FILE)
        build_cmd = BLAST_BUILD_CMD_FMT.format(
            work_dir=temp_dir, sequence_file=SEQUENCE_FILE, db_file=DB_FILE
        )
        complete = subprocess.run(
            build_cmd.split(),
            check=True,
            timeout=BLAST_TIMEOUT,
            universal_newlines=True,
        )
        complete.check_returncode()
        for db_file in temp_dir.glob(f"{DB_FILE}*"):
            _LOGGER.debug(f"Copying {db_file} to {blast_dir}.")
            shutil.copy(db_file.resolve(), blast_dir.resolve())


def run_query(sequence_file, blast_dir):
    """Build the BLAST database.

    :param str sequence_file:  path to FASTA file with sequences for DB
    :param str blast_dir:  directory for BLAST results
    """
    # The temporary directory fixes problems with spaces in paths
    with TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        _LOGGER.debug(f"Using temporary directory {temp_dir} to run BLAST.")
        shutil.copy(sequence_file.resolve(), Path(temp_dir) / SEQUENCE_FILE)
        for db_file in blast_dir.glob(f"{DB_FILE}*"):
            shutil.copy(db_file.resolve(), temp_dir.resolve())
        query_cmd = BLAST_QUERY_CMD_FMT.format(
            work_dir=temp_dir,
            sequence_file=SEQUENCE_FILE,
            db_file=DB_FILE,
            output_file=OUTPUT_FILE,
        )
        complete = subprocess.run(
            query_cmd.split(),
            check=True,
            timeout=BLAST_TIMEOUT,
            universal_newlines=True,
        )
        complete.check_returncode()
        output_file = Path(temp_dir) / Path(OUTPUT_FILE)
        _LOGGER.debug(f"Copying {output_file} to {blast_dir}.")
        shutil.copy(output_file.resolve(), blast_dir.resolve())


def process_blast(
    blast_dir, identity_cutoff, similarity_cutoff
) -> pd.DataFrame:
    """Process BLAST results.

    :param str blast_dir:  directory with BLAST results.
    :param float identity_cutoff:  keep hit if either identity or similarity
        are above cutoff
    :param float similarity_cutoff:  keep hit if either identity or similarity
        are above cutoff
    :returns:  DataFrame of results
    """
    rows = []
    with open(Path(blast_dir) / Path("output.xml")) as blast_file:
        for record in NCBIXML.parse(blast_file):
            query_id = record.query.split()[0]
            query_length = record.query_length
            for alignment in record.alignments:
                words = alignment.title.split("|")
                pdb_id = words[1]
                chain_id = words[2].split()[0]
                align_id = f"{pdb_id}_{chain_id}"
                if len(alignment.hsps) > 1:
                    _LOGGER.warning(
                        f"Ignoring extra alignments for {query_id}."
                    )
                hsp = alignment.hsps[0]
                identity = hsp.identities / max(query_length, alignment.length)
                similarity = hsp.positives / max(
                    query_length, alignment.length
                )
                if (identity >= identity_cutoff) or (
                    similarity >= similarity_cutoff
                ):
                    row = {
                        "Source": query_id,
                        "Target": align_id,
                        "Query length": query_length,
                        "Alignment length": alignment.length,
                        "Alignment score": hsp.score,
                        "Alignment bits": hsp.bits,
                        "Alignment e-value": hsp.expect,
                        "Alignment identities": hsp.identities,
                        "Alignment positives": hsp.positives,
                        "Alignment gaps": hsp.gaps,
                        "Alignment frac. identity": identity,
                        "Alignment frac. similarity": similarity,
                    }
                    rows.append(row)
    df = pd.DataFrame(rows)
    df = df[df["Source"] != df["Target"]]
    df = df.reset_index(drop=True)
    return df


def run_blast(
    fasta_path, blast_dir, identity_cutoff, similarity_cutoff, save_output=None
) -> pd.DataFrame:
    """Run all-vs-all BLAST on FASTA files.

    :param str fasta_path:  path to FASTA file with sequences
    :param str blast_dir:  directory with BLAST database
    :param float identity_cutoff:  keep hit if either identity or similarity
        are above cutoff
    :param float similarity_cutoff:  keep hit if either identity or similarity
        are above cutoff
    :param Path save_output:  path to save XML-format BLAST output (can be 
        None)
    :returns:  DataFrame with BLAST results
    """
    fasta_path = Path(fasta_path)
    blast_dir = Path(blast_dir)
    _LOGGER.info(
        f"Querying BLAST database in {blast_dir} with sequences from "
        f"{fasta_path}."
    )
    run_query(fasta_path, blast_dir)
    if save_output is not None:
        output_path = blast_dir / Path("output.xml")
        save_output = Path(save_output)
        _LOGGER.info(f"Saving raw BLAST output to {save_output}.")
        shutil.copy(output_path.resolve(), save_output.resolve())
    _LOGGER.info(f"Processing BLAST results in {blast_dir}.")
    df = process_blast(blast_dir, identity_cutoff, similarity_cutoff)
    _LOGGER.info(f"Found {len(df)} matches.")
    return df
