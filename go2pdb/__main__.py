"""Process PDB entries with specific Gene Ontology annotations."""
import logging
import argparse
from pathlib import Path
import pandas as pd
import pdbx
from . import goa, pdb, blast
from . import (
    CIF_DIR,
    GOA_DIR,
    PDB_DIR,
    FASTA_DIR,
    BLAST_DIR,
)


_LOGGER = logging.getLogger(__name__)
GO_IDS = ["GO:0090729", "GO:0031640"]
BLAST_CUTOFF = 0.9
OUTPUT_FILE = "results.xlsx"


def build_parser() -> argparse.ArgumentParser:
    """Build an argument parser."""
    parser = argparse.ArgumentParser(
        "Fetch and process a PDB GOA file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    local_group = parser.add_argument_group("local", "Local storage options")
    local_group.add_argument(
        "--blast-dir",
        help="Directory for BLAST database and results",
        default=BLAST_DIR,
    )
    local_group.add_argument(
        "--cif-dir",
        help="Directory for CIF-format file downloads",
        default=CIF_DIR,
    )
    local_group.add_argument(
        "--fasta-dir",
        help="Directory for FASTA-format file downloads",
        default=FASTA_DIR,
    )
    local_group.add_argument(
        "--goa-dir",
        help="Path to directory with GZIPped downloaded GOA file and summary of GO hits",
        default=GOA_DIR,
    )
    local_group.add_argument(
        "--pdb-dir",
        help="Directory for PDB-format file downloads",
        default=PDB_DIR,
    )
    flow_group = parser.add_argument_group(
        "flow", "Options to change workflow"
    )
    flow_group.add_argument(
        "--skip-goa",
        help="Skip download of the GOA file",
        action="store_true",
    )
    flow_group.add_argument(
        "--clobber-goa",
        help=(
            "Download the GOA file even if it is not newer than the "
            "current file"
        ),
        action="store_true",
    )
    flow_group.add_argument(
        "--skip-extract",
        help="Skip extract/download of PDB entries identified from GOA",
        action="store_true",
    )
    flow_group.add_argument(
        "--skip-blast",
        help="Skip all vs. all BLAST analysis of PDB sequences",
        action="store_true",
    )
    flow_group.add_argument(
        "--clobber-blast-db",
        help="Rebuild BLAST database even if no new sequences",
        action="store_true",
    )
    flow_group.add_argument(
        "--clobber-blast-results",
        help="Recalculate BLAST results even if no new sequences",
        action="store_true",
    )
    flow_group.add_argument(
        "--clobber-blast-analysis",
        help="Reprocess BLAST results even if no new sequences",
        action="store_true",
    )
    parser.add_argument(
        "--go-codes",
        help="GO codes for searching the PDB",
        default=GO_IDS,
        nargs="+",
    )
    parser.add_argument(
        "--identity-cutoff",
        help="cutoff for identity score (calculated via Jaccard index)",
        default=BLAST_CUTOFF,
        type=float,
    )
    parser.add_argument(
        "--output-file",
        help="file for Excel-format output",
        default=OUTPUT_FILE
    )
    parser.add_argument(
        "--log-level",
        help="Output verbosity",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
    )
    return parser


def enrich(clusters, identity_cutoff, cif_dir, goa_dir) -> pd.DataFrame:
    """Add information about PDB sequences to DataFrame.

    :param pd.DataFrame clusters:  dataframe with cluster information
    :param float identity_cutoff:  cutoff value used for identity-based clustering
    :param str cif_dir:  path to directory with CIF data
    :param str cif_dir:  path to directory with GOA data
    :returns:  dataframe with clusters enriched with information about entries
    """
    rows = []
    cif_path = Path(cif_dir) / Path(pdb.SUMMARY_FILE)
    cif_df = pd.read_excel(cif_path)
    goa_path = Path(goa_dir) / Path(goa.SUMMARY_FILE)
    goa_df = pd.read_excel(goa_path)
    goa_df["DB_Object_ID"] = goa_df["DB_Object_ID"].str.upper()
    pdb_df = goa_df["DB_Object_ID"].str.split("_", expand=True)
    goa_df["PDB ID"] = pdb_df[0]
    goa_df["Chain ID"] = pdb_df[1]
    df = pd.merge(goa_df, clusters, on="PDB ID", how="left")
    df = pd.merge(df, cif_df, on="PDB ID", how="left")
    df = df[[
        'Cluster',
        'Dep date',
        'PDB ID',
        'Chain ID',
        'Description',
        'Title',
        'Organism',
        'Experimental method',
        'Resolution (A)',
        'Date',
        'Qualifiers',
        'GO Label',
        'GO Identifier',
        'DB:Reference',
        'Evidence',
        'With',
        'Aspect',
        'Taxon_ID',
        'Assigned_By']]
    df = df.rename({
        "Cluster": f"{identity_cutoff}-identity cluster",
        'Dep date': "PDB deposit date",
        'Description': 'PDB description',
        'Title': 'PDB title',
        'Organism': 'PDB organism',
        'Experimental method': 'PDB experimental method',
        'Resolution (A)': 'PDB structure resolution (A)',
        'Date': 'GO annotation date',
        'Qualifiers': 'GO annotation qualifiers',
        'GO Label': 'GO label',
        'GO Identifier': 'GO identifier',
        'DB:Reference': 'GO database reference',
        'Evidence': 'GO annotation evidence',
        'With': 'GO evidence codes',
        'Aspect': 'GO aspect',
        'Taxon_ID': 'GO species taxon ID',
        'Assigned_By': 'GO annotation assigned_by'
    }, axis=1)
    return df


def main():
    """Main driver."""
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level, "INFO"))
    if not args.skip_goa:
        _LOGGER.info("Checking for GOA database updates.")
        goa.goa_check_fetch(
            local_dir=args.goa_dir, clobber=args.clobber_goa,
        )
    if not args.skip_extract:
        _LOGGER.info(f"Searching GOA database file for {args.go_codes}.")
        df = goa.extract(local_dir=args.goa_dir, go_codes=args.go_codes)
        pdb_ids = set(df["DB_Object_ID"])
        _LOGGER.info(f"Found {len(pdb_ids)} results.")
        _LOGGER.info(f"Checking/fetching PDB-format files in {args.pdb_dir}.")
        pdb.pdb_check_fetch(pdb_ids, args.pdb_dir)
        pdb.cif_check_fetch(pdb_ids, args.cif_dir)
        pdb.fasta_check_fetch(pdb_ids, args.fasta_dir)
        pdb.summarize_cif(pdb_ids, args.cif_dir)
    if not args.skip_blast:
        clusters = blast.run_blast(
            args.fasta_dir,
            args.blast_dir,
            args.identity_cutoff,
            clobber_db=args.clobber_blast_db,
            clobber_blast=args.clobber_blast_results,
            clobber_analysis=args.clobber_blast_analysis,
        )
        _LOGGER.info("Adding structural information to clusters.")
        clusters = enrich(clusters, args.identity_cutoff, args.cif_dir, args.goa_dir)
        _LOGGER.info(f"Writing results to {args.output_file}.")
        clusters.to_excel(args.output_file, index=False)


if __name__ == "__main__":
    main()
