import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity, get_order_residues, Silhouette, find_optimal_score_P, jaccard_index
import numpy as np

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

print(sys.argv)
active_sites = read_active_sites(sys.argv[2])
get_order_residues(active_sites)
# print one active site
# print(active_sites[0].name)

# Get residues of the active site
# Print out the first residue type of the first active site
# print((active_sites[4].residues, active_sites[5].residues))
# print((active_sites[4].newresidues, active_sites[5].newresidues))


# print('testing scoring matrix')

# print(compute_similarity(active_sites[0], active_sites[1]))



# Test clustering 1,2,3
# cluster_by_partitioning([active_sites[0], active_sites[1], active_sites[2]])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites, 2)
    write_clustering(sys.argv[3], clustering[0])

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites, 2)
    write_mult_clusterings(sys.argv[3], clusterings[0])

# Silhouette(cluster_by_partitioning(active_sites, 2)[1])
# find_optimal_score_P(active_sites)

jaccard_index(cluster_by_partitioning(active_sites, 2)[0],cluster_hierarchically(active_sites, 2)[0])

