"""Search UniProt for GO codes."""
import logging
import requests
import pandas as pd


_LOGGER = logging.getLogger(__name__)
SEARCH_GO_URL = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name&format=tsv&query=pdb%20go%3A{go}"
SEARCH_ID_URL = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name&format=tsv&query=%28accession%3A{id}%29"
MAPPING_URL = "https://rest.uniprot.org/idmapping/run"
MAPPING_DATA_URL = "https://rest.uniprot.org/idmapping/stream/{jobId}?format=tsv"


def search_id(uniprot_ids, ssl_verify) -> pd.DataFrame:
    """Search UniProt by ID.

    :param list uniprot_ids:  list of UniProt IDs to search
    :param bool ssl_verify:  does SSL work?
    :returns:  UniProt IDs and additional information
    """
    rows = []
    for uni_id in uniprot_ids:
        search_url = SEARCH_ID_URL.format(id=uni_id)
        req = requests.get(search_url, verify=ssl_verify)
        lines = req.text.splitlines()
        rows += [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(
        data=rows,
        columns=[
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt protein names",
        ],
    )
    return df.drop_duplicates(ignore_index=True)


def search_go(go_codes, ssl_verify) -> pd.DataFrame:
    """Search UniProt by GO code.

    :param list go_codes:  list of GO codes to search
    :param bool ssl_verify:  does SSL work?
    :returns:  UniProt IDs and additional information
    """
    rows = []
    for go_code in go_codes:
        _, code = go_code.split(":")
        search_url = SEARCH_GO_URL.format(go=code)
        req = requests.get(search_url, verify=ssl_verify)
        lines = req.text.splitlines()
        rows += [line.split("\t") + [go_code] for line in lines[1:]]
    df = pd.DataFrame(
        data=rows,
        columns=[
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt protein names",
            "UniProt GO code",
        ],
    )
    return df.drop_duplicates(ignore_index=True)


def get_pdb_ids(uniprot_ids, ssl_verify) -> pd.DataFrame:
    """Get PDB IDs corresponding to UniProt IDs.

    :param list uniprot_ids:  list of UniProt IDs.
    :param bool ssl_verify:  does SSL work?
    :returns:  mapping of UniProt IDs to PDB IDs.
    """
    
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "PDB",
        "ids": " ".join(uniprot_ids),
    }
    
    req = requests.post(MAPPING_URL, data=params, verify=ssl_verify)
    req_json = req.json()
    job_id = req_json["jobId"]
    
    data_url = MAPPING_DATA_URL.format(jobId=job_id)
    data_req = requests.get(data_url)
    
    lines = data_req.text.splitlines()
    rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(data=rows, columns=["UniProt entry ID", "PDB ID"])
    return df.drop_duplicates(ignore_index=True)
