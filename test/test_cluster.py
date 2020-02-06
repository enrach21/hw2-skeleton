from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_minval():
    matrix=[[1, 2, 3],[4, 5, 7],[7, 8, 10]]
    distance_matrix = [[0]*len(matrix) for i in  range(len(matrix))]
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            distance_matrix[i][j] = Euclidean_distance(matrix[i],matrix[j])
    column=0
    visited=[0]
    assert((cluster.min_val(dist, column, visited))==(5.830951894845301, 0, 1))
    column=1
    visited=[]
    assert((cluster.min_val(dist, column, visited))==(0, 1, 1))
    
def test_get_order_residues():
    filename_a = os.path.join("data", "46495.pdb")
    activesite_a = io.read_active_site(filename_a)
    assert (activesite_a.newresidues) == []
    cluster.get_order_residues(activesite_a)
    assert (activesite_a.newresidues) == ['ASP' 'LYS' 'SER' 'ARG' 'ASP' 'ASP' 'ASP']


def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert cluster.compute_similarity(activesite_a, activesite_b) == 0.0

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []
