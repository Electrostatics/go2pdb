"""Manage BLAST calculations on sequences."""
import logging
import subprocess
import shutil
import json
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
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
CLUSTER_FILE = "clusters.json"
MEMBER_FILE = "members.xlsx"


def read_fasta_dir(fasta_dir) -> list:
    """Read all of the FASTA files in a directory.

    :param str fasta_dir:  path to directory with ``*.fasta`` files.
    :returns:  list of :class:`Bio.SeqRecord` sequences
    """
    _LOGGER.debug(f"Reading FASTA records from {fasta_dir}.")
    fasta_dir = Path(fasta_dir)
    records = {}
    for fasta_file in fasta_dir.glob("*.fasta"):
        fasta_records = SeqIO.parse(fasta_file, "fasta")
        for record in fasta_records:
            if record.description not in records:
                records[record.description] = record
    _LOGGER.debug(f"Found {len(records)} unique records.")
    return list(records.values())


def new_fasta(fasta_dir, blast_dir) -> bool:
    """Check whether the FASTA results are newer than the BLAST database.

    :param str fasta_dir:  path to directory with ``*.fasta`` files.
    :param str blast_dir:  directory for BLAST results
    :returns:  ``True`` if there are new FASTA results
    """
    fasta_dir = Path(fasta_dir)
    fasta_mtime = max([f.stat().st_mtime for f in fasta_dir.glob("*.fasta")])
    blast_dir = Path(blast_dir)
    blast_mtimes = [f.stat().st_mtime for f in blast_dir.glob(f"{DB_FILE}*")]
    if not blast_mtimes:
        _LOGGER.debug("No BLAST DB files present.")
        return True
    blast_mtime = max(blast_mtimes)
    _LOGGER.debug(f"FASTA mtime = {fasta_mtime}.")
    _LOGGER.debug(f"BLAST mtime = {blast_mtime}.")
    return fasta_mtime > blast_mtime


def build_db(sequence_file, blast_dir):
    """Build the BLAST database.

    :param str sequence_file:  path to FASTA file with sequences for DB
    :param str blast_dir:  directory for BLAST results
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


def process_blast(blast_dir) -> pd.DataFrame:
    """Process BLAST results.

    :param str blast_dir:  directory with BLAST results.
    :returns:  DataFrame of results
    """
    results = {}
    with open(Path(blast_dir) / Path("output.xml")) as blast_file:
        for record in NCBIXML.parse(blast_file):
            query_id = record.query.split("|")
            pdb_id, chain_num = query_id[0].split("_")
            query_id = f"pdb|{pdb_id}|{chain_num}|" + "|".join(query_id[1:])
            query_length = record.query_length
            query_dict = results.get(query_id, {})
            for alignment in record.alignments:
                align_id = alignment.title
                if len(alignment.hsps) > 1:
                    _LOGGER.warning(
                        f"Ignoring extra alignments for {query_id}."
                    )
                hsp = alignment.hsps[0]
                score = hsp.identities / max(query_length, alignment.length)
                query_dict[align_id] = score
            _LOGGER.debug(f"Found {len(query_dict)} matches for {query_id}.")
            results[query_id] = query_dict
    return pd.DataFrame(results)


def cluster_index(value, clusters) -> int:
    """Check whether something is in a cluster.

    :param value:  value to check
    :param list clusters:  list of clusters (set)
    :returns:  index of cluster containing the value
    :raises IndexError:  if value not in any cluster.
    """
    for icluster, cluster in enumerate(clusters):
        if value in cluster:
            return icluster
    raise IndexError(value)


def cluster(df) -> list:
    """Cluster DataFrame by grouping non-NA entries into clusters.

    :param pd.DataFrame df:  DataFrame to cluster
    :returns:  list of clusters (sets)
    """
    clusters = []
    for hit1, row in df.iterrows():
        try:
            icluster = cluster_index(hit1, clusters)
        except IndexError:
            clusters.append({hit1})
            icluster = len(clusters) - 1
        row = row.dropna()
        for hit2 in row.index:
            clusters[icluster].add(hit2)
    return clusters


def transform_clusters(cluster_list) -> tuple:
    """Transform clusters into a membership dictionary and component dictionary.

    :param list cluster_list:  list of clusters (sets)
    :returns: tuple with (component dict, membership dict)
    """
    cluster_components = {}
    cluster_membership = {}
    for clust in cluster_list:
        clust = sorted(list(clust))
        clust_name = clust[0].split("|")[4]
        count_name = list(cluster_components.keys()).count(clust_name)
        if count_name > 1:
            count_name = f"{clust_name} ({count_name})"
        for hit in clust:
            members = cluster_membership.get(hit, set())
            members.add(clust_name)
            cluster_membership[hit] = members
        cluster_components[clust_name] = clust
    _LOGGER.info(f"Found {len(cluster_list)} clusters.")
    for seq, clust in cluster_membership.items():
        if len(clust) > 1:
            _LOGGER.warning(f"{seq} appears in {len(clust)} clusters: ")
            for name in clust:
                _LOGGER.warning(f"    {name}")
    return cluster_components, cluster_membership


def run_blast(
    fasta_dir,
    blast_dir,
    identity_cutoff,
    clobber_db=False,
    clobber_blast=False,
    clobber_analysis=False,
) -> pd.DataFrame:
    """Run all-vs-all BLAST on FASTA files.

    :param str fasta_dir:  directory with FASTA ``*.fasta`` files
    :param str blast_dir:  directory for BLAST results
    :param float identity_cutoff:  cutoff for accepting BLAST results based on Jaccard identity
    :param bool clobber_db:  rebuild database
    :param bool clobber_blast:  re-run BLAST
    :param bool clobbber_analysis:  re-process BLAST data
    :returns:  DataFrame assigning sequences to clusters
    """
    fasta_dir = Path(fasta_dir)
    blast_dir = Path(blast_dir)
    new_input = new_fasta(fasta_dir, blast_dir)
    _LOGGER.debug(f"New input? {new_input}")
    _LOGGER.debug(f"Clobber DB? {clobber_db}")
    _LOGGER.debug(f"Clobber BLAST? {clobber_blast}")
    _LOGGER.debug(f"Clobber processed results? {clobber_analysis}")
    if new_input or clobber_db:
        _LOGGER.info(
            f"Building BLAST database in {blast_dir} with sequences from {fasta_dir}."
        )
        sequences = read_fasta_dir(fasta_dir)
        sequence_file = blast_dir / Path(SEQUENCE_FILE)
        _LOGGER.debug(f"Writing sequence data to {sequence_file}.")
        SeqIO.write(sequences, sequence_file, "fasta")
        build_db(sequence_file, blast_dir)
    if new_input or clobber_db or clobber_blast:
        _LOGGER.info(
            f"Querying BLAST database in {blast_dir} with sequences from {fasta_dir}."
        )
        run_query(sequence_file, blast_dir)
    if new_input or clobber_blast or clobber_blast or clobber_analysis:
        _LOGGER.info(f"Processing BLAST results in {blast_dir}.")
        df = process_blast(blast_dir)
        tot_match = df.shape[0] * df.shape[1]
        _LOGGER.debug(
            f"Found {df.count().sum()} significant matches (of {tot_match} possible)."
        )
        df[df < identity_cutoff] = np.nan
        _LOGGER.info(
            f"Found {df.count().sum()} matches above cutoff ({identity_cutoff})."
        )
        cluster_list = cluster(df)
        cluster_components, cluster_membership = transform_clusters(
            cluster_list
        )
        cluster_path = blast_dir / Path(CLUSTER_FILE)
        _LOGGER.info(f"Writing clusters to {cluster_path}.")
        with open(cluster_path, "wt") as cluster_file:
            json.dump(cluster_components, cluster_file, indent=2)
        rows = []
        for hit, clusters in cluster_membership.items():
            for clust in clusters:
                words = hit.split("|")
                pdb_id = words[1]
                chain_num = int(words[2])
                chains = words[3]
                description = words[4]
                organism = words[5]
                rows.append(
                    (clust, pdb_id, chain_num, chains, description, organism)
                )
        member_df = pd.DataFrame(
            rows,
            columns=[
                "Cluster",
                "PDB ID",
                "Chain #",
                "Chains",
                "Description",
                "Organism",
            ],
        )
        member_path = blast_dir / Path(MEMBER_FILE)
        _LOGGER.info(f"Writing sequence cluster membership to {member_path}.")
        member_df.to_excel(member_path, index=False)
    return member_df
