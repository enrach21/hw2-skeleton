from hw2skeleton import cluster
from hw2skeleton import io
import pandas as pd
import os
import glob
import numpy as np


# Make a dictionary for the various combination of amino acid changes
AminoAcid = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE','PRO', 'SER', 'THR', 'TRP','TYR','VAL','ASX','XLE', 'GLX','XAA', '*']
BLOSUM = glob.glob('BLOSUM80.txt')
data = pd.read_table(BLOSUM[0], sep=" ", header = None, names = AminoAcid)
data = data.set_index([pd.Index(AminoAcid)])
print(data)
Dict = {}
for i in AminoAcid:
    for j in AminoAcid:
        Dict[(i,j)]=data[i][j]

def test_minval():
    matrix=[[1, 2, 3],[4, 5, 7],[7, 8, 10]]
    distance_matrix = [[0]*len(matrix) for i in  range(len(matrix))]
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            distance_matrix[i][j] = cluster.Euclidean_distance(matrix[i],matrix[j])
    column=0
    visited=[0]
    assert ((cluster.min_val(distance_matrix, column, visited)) ==(5.830951894845301, 0, 1))
    column=1
    visited=[]
    assert ((cluster.min_val(distance_matrix, column, visited)) ==(0, 1, 1))

def test_get_order_residues():
    filename_a = os.path.join("data", "46495.pdb")
    activesite_a = io.read_active_site(filename_a)
    assert np.array_equal((activesite_a.newresidues),[])
    cluster.get_order_residues([activesite_a])
    assert np.array_equal((activesite_a.newresidues),['ASP', 'LYS', 'SER', 'ARG', 'ASP', 'ASP', 'ASP'])



def test_similarity():
    filename_a = os.path.join("data", "46495.pdb")
    filename_b = os.path.join("data", "23812.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    cluster.get_order_residues([activesite_a])
    cluster.get_order_residues([activesite_b])
    # update this assertion
    assert (cluster.compute_similarity(activesite_a, activesite_b)) == 13

def test_dist():
    filename_a = os.path.join("data", "46495.pdb")
    filename_b = os.path.join("data", "23812.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    cluster.get_order_residues([activesite_a])
    cluster.get_order_residues([activesite_b])
    sim1=cluster.compute_similarity(activesite_a, activesite_a)
    sim2=cluster.compute_similarity(activesite_a, activesite_a)
    sim3=cluster.compute_similarity(activesite_a, activesite_b)
    sim4=cluster.compute_similarity(activesite_b, activesite_a)
    assert cluster.Euclidean_distance(sim1, sim2) == 0
    assert cluster.Euclidean_distance(sim3, sim4)== cluster.Euclidean_distance(sim4, sim3)
    assert cluster.Euclidean_distance(sim3, sim4) > 0



def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701, 10701,10814,13052,14181,15813]


    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert len(cluster.cluster_by_partitioning(active_sites,2)[0]) == 2

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701, 10701,10814,13052,14181,15813]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert len(cluster.cluster_hierarchically(active_sites,2)[0]) == 2




