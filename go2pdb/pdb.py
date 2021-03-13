"""Routines for fetching things from the PDB website."""
import logging
import json
from os import access
from pathlib import Path
from datetime import date
import requests
import pandas as pd



_LOGGER = logging.getLogger(__name__)
SUMMARY_FILE = "summary.xlsx"
PDB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
PDB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v1/query?json={query}"
CHUNK_SIZE = 180
GRAPHQL_URL = "https://data.rcsb.org/graphql"
GRAPHQL_QUERY = (
    "{{"
    "  entries(entry_ids: {pdb_ids}) {{"
    "    rcsb_id"
    "    rcsb_accession_info {{ deposit_date }}"
    "    struct {{ title pdbx_descriptor }}"
    "    exptl {{ method }}"
    "    refine {{ ls_d_res_high }}"
    "    polymer_entities {{"
    "      rcsb_id"
    "      entity_poly {{"
    "        pdbx_strand_id"
    "        rcsb_entity_polymer_type"
    "        pdbx_seq_one_letter_code_can"
    "      }}"
    "      uniprots {{ rcsb_id }}"
    "    }}"
    "  }}"
    "}}"
)



def keyword_search(keyword, ssl_verify) -> pd.DataFrame:
    """Perform a keyword search of the PDB.

    :param str keyword:  keyword to search
    :param bool ssl_verify:  does SSL work?
    :returns:  matching PDB IDs
    """
    pdb_ids = set()
    for key in (keyword, keyword.upper(), keyword.lower()):
        for field in ("struct_keywords.text", "struct_keywords.pdbx_keywords"):
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": field,
                        "operator": "contains_words",
                        "value": key
                    }
                },
                "request_options": {
                    "return_all_hits": True
                },
                "return_type": "entry"
            }
            search_url = PDB_SEARCH_URL.format(query=json.dumps(query))
            req = requests.get(search_url, verify=ssl_verify)
            req_dict = req.json()
            for result in req_dict["result_set"]:
                pdb_id = result["identifier"].upper()
                pdb_ids.add(pdb_id)
    df = pd.DataFrame(data=list(pdb_ids), columns=["PDB ID"])
    df["PDB keyword"] = keyword
    return df


def metadata(pdb_ids, ssl_verify) -> pd.DataFrame:
    """Get metadata for PDB IDs.

    :param list pdb_ids:  list of PDB IDs
    :param bool ssl_verify:  does SSL work?
    :returns:  DataFrame with metadata
    """
    pdb_ids = sorted(list(pdb_ids))
    rows = []
    for id_list in [pdb_ids[i:i+CHUNK_SIZE] for i in range(0, len(pdb_ids), CHUNK_SIZE)]:
        _LOGGER.debug(f"Fetching metadata for {id_list}.")
        query = GRAPHQL_QUERY.format(pdb_ids=id_list)
        query = query.replace("'", '"')
        req = requests.get(GRAPHQL_URL, params={"query": query}, verify=ssl_verify)
        results = req.json()
        for result in results["data"]["entries"]:
            pdb_id = result.pop("rcsb_id")
            struct = result.pop("struct")
            descriptor = struct.pop("pdbx_descriptor")
            accession = result.pop("rcsb_accession_info")
            dep_date = accession.pop("deposit_date")
            dep_date, _ = dep_date.split("T")
            dep_date = date.fromisoformat(dep_date)
            title = struct.pop("title")
            experiments = result.pop("exptl")
            if len(experiments) > 1:
                _LOGGER.warning(f"Only using first experiment of {pdb_id} for annotation.")
            experiment = experiments[0]
            method = experiment.pop("method")
            refinements = result.pop("refine")
            if refinements is not None:
                if len(refinements) > 1:
                    _LOGGER.warning(f"Only using first refinement of {pdb_id} for annotation.")
                refinement = refinements[0]
                resolution = refinement.pop("ls_d_res_high")
            else:
                resolution = None
            for polymer in result.pop("polymer_entities"):
                chain_id = polymer.pop("rcsb_id")
                entity = polymer.pop("entity_poly")
                strand_ids = entity.pop("pdbx_strand_id")
                strand_type = entity.pop("rcsb_entity_polymer_type")
                sequence = entity.pop("pdbx_seq_one_letter_code_can")
                uniprots = polymer.pop("uniprots")
                if uniprots is not None:
                    if len(uniprots) > 1:
                        _LOGGER.warning(f"Only using first UniProt ID of {pdb_id} for annotation")
                    uniprot = uniprots[0].pop("rcsb_id")
                else:
                    uniprot = None
                row = {"PDB ID": pdb_id, "PDB deposit date": dep_date, "PDB method": method, "PDB resolution (A)": resolution, "PDB description": descriptor, "PDB title": title, "PDB chain ID": chain_id, "PDB strand ID(s)": strand_ids, "PDB strand type": strand_type, "PDB strand sequence": sequence, "PDB strand UniProt": uniprot}
                rows.append(row)
    return pd.DataFrame(rows)