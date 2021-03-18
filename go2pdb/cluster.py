"""Cluster sequences/structures based on BLAST results."""
import logging
from hashlib import sha1
from subprocess import NORMAL_PRIORITY_CLASS
import pandas as pd


_LOGGER = logging.getLogger(__name__)


CLUSTER_FILE = "clusters.json"
MEMBER_FILE = "members.xlsx"


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
    :returns:  list of clusters (lists)
    """
    clusters = []
    for _, row in df.iterrows():
        source = row["Source"]
        target = row["Target"]
        try:
            icluster = cluster_index(source, clusters)
        except IndexError:
            clusters.append({source})
            icluster = len(clusters) - 1
        clusters[icluster].add(target)
    clusters = [sorted(list(c)) for c in clusters]
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def transform_clusters(cluster_list) -> pd.DataFrame:
    """Transform clusters into a membership dictionary and component dictionary.

    :param list cluster_list:  list of clusters (sets)
    :returns: DataFrame with chains mapped to clusters
    """
    components = {}
    chains = set()
    for cluster in cluster_list:
        cluster = sorted(cluster)
        chains |= set(cluster)
        key = cluster[0]
        components[key] = cluster
    rows = []
    for key in components:
        for chain in sorted(components[key]):
            rows.append((key, chain))
    return pd.DataFrame(rows, columns=["Cluster", "Chain"])

    # tot_match = df.shape[0] * df.shape[1]
    # _LOGGER.debug(
    #     f"Found {df.count().sum()} significant matches (of {tot_match} possible)."
    # )
    # df[df < identity_cutoff] = np.nan
    # _LOGGER.info(
    #     f"Found {df.count().sum()} matches above cutoff ({identity_cutoff})."
    # )
    # cluster_list = cluster(df)
    # cluster_components, cluster_membership = transform_clusters(
    #     cluster_list
    # )
    # cluster_path = blast_dir / Path(CLUSTER_FILE)
    # _LOGGER.info(f"Writing clusters to {cluster_path}.")
    # with open(cluster_path, "wt") as cluster_file:
    #     json.dump(cluster_components, cluster_file, indent=2)
    # rows = []
    # for hit, clusters in cluster_membership.items():
    #     for clust in clusters:
    #         words = hit.split("|")
    #         pdb_id = words[1]
    #         chain_num = int(words[2])
    #         chains = words[3]
    #         description = words[4]
    #         organism = words[5]
    #         rows.append(
    #             (clust, pdb_id, chain_num, chains, description, organism)
    #         )
    # member_df = pd.DataFrame(
    #     rows,
    #     columns=[
    #         "Cluster",
    #         "PDB ID",
    #         "Chain #",
    #         "Chains",
    #         "Description",
    #         "Organism",
    #     ],
    # )
    # member_path = blast_dir / Path(MEMBER_FILE)
    # _LOGGER.info(f"Writing sequence cluster membership to {member_path}.")
    # member_df.to_excel(member_path, index=False)
    # return member_df
