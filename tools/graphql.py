"""Test GraphQL."""
import requests
from pprint import pprint


GRAPHQL_URL = "https://data.rcsb.org/graphql"
GRAPHQL_QUERY = (
    "{{"
    "  entries(entry_ids: {pdb_ids}) {{"
    "    rcsb_id"
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


def main():
    """Main driver."""
    # query = '{ entries(entry_ids: ["4HHB", "12CA", "3PQR"]) { exptl { method } } }'
    query = GRAPHQL_QUERY.format(
        pdb_ids=["1FAS", "7BLO", "1MAH", "2AAI", "2OS6"]
    )
    query = query.replace("'", '"')
    params = {"query": query}
    req = requests.get(GRAPHQL_URL, params=params)
    results = req.json()
    for result in results["data"]["entries"]:
        pdb_id = result.pop("rcsb_id")
        struct = result.pop("struct")
        descriptor = struct.pop("pdbx_descriptor")
        title = struct.pop("title")
        experiments = result.pop("exptl")
        if len(experiments) > 1:
            print("Only using first experiment for annotation.")
        experiment = experiments[0]
        method = experiment.pop("method")
        refinements = result.pop("refine")
        if refinements is not None:
            if len(refinements) > 1:
                print("Only using first refinement for annotation.")
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
                    print("Only using first UniProt ID for annotation")
                uniprot = uniprots[0].pop("rcsb_id")
            else:
                uniprot = None
            row = {
                "PDB ID": pdb_id,
                "PDB method": method,
                "PDB resolution (A)": resolution,
                "PDB description": descriptor,
                "PDB title": title,
                "PDB chain ID": chain_id,
                "PDB strand ID(s)": strand_ids,
                "PDB strand type": strand_type,
                "PDB strand sequence": sequence,
                "PDB strand UniProt": uniprot,
            }
            pprint(row)
        print(result)


if __name__ == "__main__":
    main()
