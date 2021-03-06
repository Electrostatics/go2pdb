"""Fetch and process PDB Gene Ontology Association data."""
import logging
import gzip
from ftplib import FTP
from pathlib import Path
from datetime import datetime, date
import pandas as pd


_LOGGER = logging.getLogger(__name__)
NOW = datetime.now()
FTP_SERVER = "ftp.ebi.ac.uk"
FTP_DIR = "pub/databases/GO/goa/PDB"
FTP_FILENAME = "goa_pdb.gaf.gz"
SUMMARY_FILE = "goa_summary.xlsx"
GOA_COMMENT = "!"
GOA_COLUMNS = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifiers",
    "GO Identifier",
    "DB:Reference",
    "Evidence",
    "With",
    "Aspect",
    "DB_Object_Name",
    "Synonym",
    "DB_Object_Type",
    "Taxon_ID",
    "Date",
    "Assigned_By",
]
PDB_COLUMN = GOA_COLUMNS.index("DB_Object_ID")
QUAL_COLUMN = GOA_COLUMNS.index("Qualifiers")
GO_COLUMN = GOA_COLUMNS.index("GO Identifier")
GOA_EVIDENCE = {
    "EXP": "inferred from experiment",
    "IMP": "inferred from mutant phenotype",
    "IC": "inferred by curator",
    "IGI": "inferred from genetic interaction",
    "IPI": "inferred from physical interaction",
    "ISS": "inferred from sequence or structural similarity",
    "IDA": "inferred from direct assay",
    "IEP": "inferred from expression pattern",
    "IEA": "inferred from electronic annotation",
    "TAS": "traceable author statement",
    "NAS": "non-traceable author statement",
    "NR": "not recorded",
    "ND": "no biological data available",
    "RCA": "inferred from reviewed computational analysis",
    "IBA": "inferred from biological aspect of ancestor",
    "ISM": "inferred from sequence model",
    "ISO": "inferred from sequence orthology",
}


def compare_mtime(ftp, gzip_file, ftp_filename) -> bool:
    """Get modification times of FTP files.

    :param FTP ftp:  connected FTP object
    :param gzip.GzipFile gzip_file:  gzip file open for reading
    :param string ftp_filename:  name of file on FTP server
    :returns:  True if ftp file is newer than local file
    """
    dir_list = []
    ftp.dir(dir_list.append)
    mtime_dict = {}
    for line in dir_list:
        words = line.split()
        month = words[5]
        day = words[6]
        if ":" in words[7]:
            year = NOW.strftime("%Y")
        else:
            year = words[7]
        dt = datetime.strptime(f"{day} {month} {year}", "%d %b %Y")
        file_date = date(dt.year, dt.month, dt.day)
        file_name = words[8]
        mtime_dict[file_name] = file_date
    gzip_file.peek(8)
    local_mtime = date.fromtimestamp(gzip_file.mtime)
    if mtime_dict[ftp_filename] > local_mtime:
        _LOGGER.debug(
            f"Remote file date {mtime_dict[ftp_filename]} is newer than "
            f"the local file date {local_mtime}."
        )
        return True
    return False


def check_fetch(
    local_path,
    ftp_server=FTP_SERVER,
    ftp_dir=FTP_DIR,
    ftp_filename=FTP_FILENAME,
):
    """Check to see if the FTP file is newer and fetch, if needed.

    :param str local_path:  directory for local gzip file
    :param str ftp_server:  FQDN of FTP server
    :param str ftp_dir:  directory of GOA file on FTP server
    :param str ftp_filename:  name of GOA file on FTP server
    """
    ftp = FTP(ftp_server)
    _LOGGER.debug(f"Logging into {ftp_dir}.")
    ftp.login()
    _LOGGER.debug(f"Changing directory to {ftp_dir}.")
    ftp.cwd(ftp_dir)
    download = None
    if local_path.exists():
        _LOGGER.debug(f"Comparing mtimes of remote and local files.")
        with gzip.GzipFile(local_path, "r") as gzip_file:
            download = compare_mtime(ftp, gzip_file, ftp_filename)
    else:
        _LOGGER.info(f"Local file {local_path} does not exist.")
        download = True
    if download:
        _LOGGER.info(f"Downloading file.")
        _LOGGER.info(f"Fetching {ftp_filename} from {ftp_server}/{ftp_dir}.")
        with open(local_path, "wb") as gzip_file:
            ftp.retrbinary(f"RETR {ftp_filename}", gzip_file.write)


def extract(local_path, go_codes) -> pd.DataFrame:
    """Extract entries with specific GO codes.

    :param str local_path:  path to local GOA gzip file
    :param list go_codes:  GO codes to search for
    :returns:  set of matching PDB IDs
    """
    _LOGGER.debug(f"Reading {local_path}.")
    rows = []
    with gzip.open(local_path, "rt") as gzip_file:
        for line in gzip_file:
            if line[0] != GOA_COMMENT:
                words = line.strip().split("\t")
                go_code = words[GO_COLUMN]
                qual = words[QUAL_COLUMN]
                if (go_code in go_codes) and ("NOT" not in qual):
                    rows.append(words)
    df = pd.DataFrame(rows, columns=GOA_COLUMNS)
    df = df.drop(
        [
            "DB",
            "DB_Object_Symbol",
            "DB_Object_Name",
            "Synonym",
            "DB_Object_Type",
        ],
        axis=1,
    )
    df["Aspect"] = df["Aspect"].replace(
        {
            "P": "biological process",
            "F": "molecular function",
            "C": "cellular component",
        }
    )
    df["Evidence"] = df["Evidence"].replace(GOA_EVIDENCE)
    df["Date"] = pd.to_datetime(df["Date"], format=r"%Y%m%d")
    df = df.rename(
        {
            "DB_Object_ID": "GOA DB object ID",
            "Qualifiers": "GOA qualifiers",
            "GO Identifier": "GOA GO code",
            "DB:Reference": "GOA DB reference",
            "Evidence": "GOA evidence",
            "With": "GOA additional evidence",
            "Aspect": "GOA aspect",
            "Taxon_ID": "GOA taxon ID",
            "Date": "GOA annotation date",
            "Assigned_By": "GOA assigned by",
        },
        axis=1,
    )
    return df
