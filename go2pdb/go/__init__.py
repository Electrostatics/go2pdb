"""Routines for handling Gene Ontology information."""
import logging
import requests


_LOGGER = logging.getLogger(__name__)
BIOENTITY_API = "http://api.geneontology.org/api/bioentity"
DB_REFS = {
    "GO_REF:0000002": "InterPro2GO",
    "GO_REF:0000043": "UniProt Keywords2GO",
    "GO_REF:0000044": "UniProt Subcellular Location2GO",
    "GO_REF:0000003": "EC2GO",
    "GO_REF:0000104": "UniRule2GO",
    "GO_REF:0000107": "Ensembl & EnsemblGenomes",
    "GO_REF:0000041": "UniPathway2GO",
    "GO_REF:0000108": "Gene Ontology Consortium",
    "GO_REF:0000115": "RNACentral",
}


def go_mapping(go_codes) -> dict:
    """Get human-readable labels for GO codes.

    :param list go_codes:  list of GO codes
    :returns:  dictionary mapping codes to labels
    """
    labels = {}
    for code in go_codes:
        url = f"{BIOENTITY_API}/{code}"
        req = requests.get(url)
        go_entry = req.json()
        labels[code] = go_entry["label"]
    return labels
