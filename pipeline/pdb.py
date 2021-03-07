"""Routines for fetching things from the PDB website."""
import logging
from pathlib import Path
from io import StringIO
from datetime import date
import requests
import pandas as pd
import pdbx
from Bio import SeqIO


_LOGGER = logging.getLogger(__name__)
SUMMARY_FILE = "summary.xlsx"
PDB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
PDB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


def pdb_check_fetch(pdb_ids, dir, url_fmt=PDB_PDB_URL):
    """Check to see if a PDB file exists and fetch it if it doesn't.

    :param list pdb_ids:  list of PDB IDs to check/fetch
    :param str dir:  path to directory where files will reside
    :param str url_fmt:  format-string for URL
    """
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        pdb_id, _ = pdb_id.split("_")
        _LOGGER.debug(f"Checking/fetching {pdb_id}.")
        path = Path(dir) / Path(f"{pdb_id}.pdb")
        if not path.exists():
            pdb_url = url_fmt.format(pdb_id=pdb_id)
            _LOGGER.info(f"Fetching {pdb_url}.")
            with requests.get(pdb_url) as req:
                with open(path, "wt") as pdb_file:
                    pdb_file.write(req.text)


def cif_check_fetch(pdb_ids, dir, url_fmt=PDB_CIF_URL):
    """Check to see if a CIF file exists and fetch it if it doesn't.

    :param list pdb_ids:  list of PDB IDs to check/fetch
    :param str dir:  path to directory where files will reside
    :param str url_fmt:  format-string for URL
    """
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        pdb_id, _ = pdb_id.split("_")
        _LOGGER.debug(f"Checking/fetching {pdb_id}.")
        path = Path(dir) / Path(f"{pdb_id}.cif")
        if not path.exists():
            pdb_url = url_fmt.format(pdb_id=pdb_id)
            _LOGGER.info(f"Fetching {pdb_url}.")
            with requests.get(pdb_url) as req:
                with open(path, "wt") as pdb_file:
                    pdb_file.write(req.text)


def fasta_check_fetch(pdb_ids, dir, url_fmt=PDB_FASTA_URL):
    """Check to see if FASTA files exist and fetch them if they don't.

    :param list pdb_ids:  list of PDB IDs to check/fetch
    :param str dir:  path to directory where files will reside
    :param str url_fmt:  format-string for URL
    """
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        pdb_id, chain_id = pdb_id.split("_")
        chain_id = chain_id.upper()
        _LOGGER.debug(f"Checking/fetching {pdb_id}.")
        path = Path(dir) / Path(f"{pdb_id}_{chain_id}.fasta")
        if not path.exists():
            pdb_url = url_fmt.format(pdb_id=pdb_id)
            _LOGGER.info(f"Fetching {pdb_url}.")
            with requests.get(pdb_url) as req:
                with open(path, "wt") as pdb_file:
                    fasta_text = req.text
            with StringIO(fasta_text) as fasta_file:
                records = list(SeqIO.parse(fasta_file, "fasta"))
            chain_record = None
            for record in records:
                words = record.description.split("|")
                chain_text = words[1]
                if "Chains" in chain_text:
                    _, chain_list = chain_text.split()
                    chain_list = [c.upper() for c in chain_list.split(",")]
                    if chain_id in chain_list:
                        chain_record = record
                        break
                elif "Chain" in chain_text:
                    _, test_chain = chain_text.split()
                    if test_chain.upper() == chain_id:
                        chain_record = record
                        break
                else:
                    err = f"Don't know how to parse description {record.description}."
                    raise ValueError(err)
            if chain_record is not None:
                SeqIO.write([chain_record], path, "fasta")
            else:
                descriptions = [rec.description for rec in records]
                err = f"Unable to find chain {chain_id} in FASTA with descriptions {descriptions}."
                raise ValueError(err)


def summarize_cif(pdb_ids, cif_dir, summary_file=SUMMARY_FILE):
    """Summarize information in CIF files for the specified PDBs and write to
    Excel file.

    :param set pdb_ids:  list of PDB IDs to summarize
    :param str cif_dir:  path to directory of CIF file
    :param str summary_file:  name of summary Excel file in CIF directory
    """
    cif_dir = Path(cif_dir)
    df_path = cif_dir / Path(summary_file)
    if df_path.exists():
        df = pd.read_excel(df_path)
        pdb_col = df["PDB ID"]
        _LOGGER.info(f"Got {df.shape} CIF summary DataFrame from {df_path}.")
    else:
        df = None
        pdb_col = pd.Series()
        _LOGGER.info("Creating new CIF summary DataFrame.")
    rows = []
    _LOGGER.info("Checking for missing CIF data.")
    pdb_ids = {pdb_id.upper().split("_")[0] for pdb_id in pdb_ids}
    for pdb_id in sorted(list(pdb_ids)):
        pdb_id = pdb_id.upper().split("_")[0]
        if pdb_id not in pdb_col:
            cif_path = cif_dir / Path(f"{pdb_id}.cif")
            with open(cif_path, "rt") as cif_file:
                container = pdbx.load(cif_file)[0]
            dep_date = container.get_object("pdbx_database_status").get_value("recvd_initial_deposition_date")
            words = dep_date.split("-")
            dep_date = date(int(words[0]), int(words[1]), int(words[2]))
            title = container.get_object("struct").get_value("title")
            method = container.get_object("exptl").get_value("method")
            try:
                resolution = container.get_object("refine").get_value("ls_d_res_high")
            except AttributeError:
                resolution = None
            row = {"PDB ID": pdb_id}
            row["Dep date"] = dep_date
            row["Title"] = title
            row["Experimental method"] = method
            row["Resolution (A)"] = resolution
            rows.append(row)
            _LOGGER.debug(f"Got CIF data for {pdb_id}.")
    new_df = pd.DataFrame(rows)
    if df is not None:
        df = df.append(new_df, ignore_index=True)
    else:
        df = new_df
    _LOGGER.info(f"Writing {df.shape} CIF summary DataFrame to {df_path}.")
    df.to_excel(df_path, index=False)
