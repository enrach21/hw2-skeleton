from .utils import Atom, Residue, ActiveSite
import numpy as np
from .io import read_active_sites, write_clustering, write_mult_clusterings, local_align
import math
from statistics import mean
import random
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import silhouette_score


def get_order_residues(active_sites):
    """
    Given a list of Active sites, add a feature 'newresidues' that is a list
    of residues 

    Input: a list of ActiveSite instances
    Output: Nothing
    """

    for active in active_sites: # a for loop going throught the active sites
        print(active.name)
        res_mean = [] # store mean coords
        res = np.array([]) # Store old res order
        res_list = np.array([]) # updated res list
        for residue in active.residues: # 
            Atom = residue.atoms
            # print(Atom)
            Dim1 = np.array([])
            Dim2 = np.array([])
            Dim3 = np.array([])
            for a in Atom:
                Dim1 = np.append(Dim1, a.coords[0])
                Dim2 = np.append(Dim2, a.coords[1])
                Dim3 = np.append(Dim3, a.coords[2])
            res_mean += [list([np.mean(Dim1), np.mean(Dim2), np.mean(Dim3)])] # Obtains the mean of all three dem
            res = np.append(res, residue.type)
        res_list = reorder(res_mean, res)
        active.newresidues = (res_list) # sets new res order to active site
    return

def reorder(res_mean, res):  
    """
    Given a list of mean coordinate values for each residue in an active site. Return
    a reordered list based on the distance between each residue. With the two residues
    closest together next to each other in the array and then add 

    Input: mean coordinate positions of residues, and list of residues before reordering
    Output: reordered residues

    """ 
    res_list = np.array([]) # updated res list         
    num_of_res = len(res_mean) # get the number of residues
    distance_matrix = [[0]*num_of_res for i in  range(num_of_res)] # Make a distance matrix of the residue positions
    for i in range(num_of_res):
        for j in range(num_of_res):
            distance_matrix[i][j] = Euclidean_distance(res_mean[i],res_mean[j])
    dist = pd.DataFrame(distance_matrix)
    # print(dist)
    visited = np.array([]) # Checks which residues have been already added to the new list
    minval = float("inf") # value used to keep track of minimum postions
    for i in range(len(dist)): # Goes throuhg distance matrix of mean residues and finds nearest ones
        col = dist[i]
        for j in range(len(dist)):
            if col[j] != 0 and minval > col[j]:
                    x = i
                    y = j
                    minval = col[j]
    dist[x][y]=float("inf") # set residues taken to inf so they won't be selected again
    dist[y][x]=float("inf")
    res_list = np.append(res_list, res[x])
    res_list = np.append(res_list, res[y])
    visited = np.append(visited , [x])
    visited  = np.append(visited , [y])
    while len(visited) < len(res_mean): # Iteratively repeat till all residues have been added to res list
        start=min_val(dist, x, visited)
        end=min_val(dist, y, visited)
        if start[0] < end[0]: # add the shortest distance next either to the begin or end of the resdiue order
            x=start[2]
            y1=start[1]
            dist[x][y1]=float("inf")
            dist[y1][x]=float("inf")
            visited  = np.insert(visited ,0, [x])
            res_list  = np.insert(res_list ,0, res[x])
        else:
            x1=end[1]
            y=end[2]
            dist[x1][y]=float("inf")
            dist[y][x1]=float("inf")
            visited  = np.append(visited , [y])
            res_list  = np.append(res_list , res[y])
    return res_list

def min_val(dist, column, visited):
    """
    Given a distance matrix find the min_value in the distance matrix and it's column positions

    Input: 
    dist is a distance matrix
    column is an index of the column of the distancer matrix being visited
    vististed is a list that keeps track of which columns have already been added to the updated list
    Output: minimum vlaue, x postion, y position

    """
    minval = float("inf")
    test = range(len(dist))
    column
    col = dist[column]
    for j in test:
        if j not in visited and minval > col[j]:
                y1 = j
                minval = col[j]
    return minval, column, y1

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances using local allignment on new residue scores

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    score = local_align(site_a.newresidues, site_b.newresidues)


    return score[0]

def Euclidean_distance(x,y):
    """
    Given too vectors calculate the euclidean distance between them

    Input: Two Vectors
    Output: Distance
    """
    D = math.sqrt(sum([(a-b) ** 2 for a, b in zip(x, y)]))
    return D

def cluster_by_partitioning(active_sites, k):
    """
    Cluster a given set of ActiveSite instances using a partitioning method (K-means).

    Input: a list of ActiveSite instances and K
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances) as well as a dataframe with info on how they were clustered
    """

    # Obtain is the updated matrix with score calculated from compute_similarity
    score_matrix = [[0]*len(active_sites) for i in  range(len(active_sites))]
    X = -1
    active = np.array([])
    for i in active_sites:
        X += 1
        Y = -1
        active = np.append(active, i.name)
        for j in active_sites:
            Y += 1
            score_matrix[X][Y] = compute_similarity(i,j) 
    # Initialization
    d = {'Scores': score_matrix, 'Cluster': [0] * (len(score_matrix)), 'Active':active}
    df = pd.DataFrame(data=d) # A dataframe of score comparison between all active sites and 
    # generate PCA to visulize all the vectors
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(score_matrix)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])
    df['principal component 1']=principalDf['principal component 1']
    df['principal component 2']=principalDf['principal component 2']
    # Select K random points based on the K
    randomK = random.sample(range(len(active_sites)), k) #Randomly select a number in the range of active sites without resampling
    centroids = [] 
    for i in randomK:
        centroids += ([score_matrix[i]]) # Get the score for the number selected and put it into a vector called centroid
    # Assign each active site to the nearest centroid
    df = assignment(df, centroids)
    # Repeat till completing the cluster
    count = 0
    while True:
        closest_centroids = df['Cluster']
        centroids = update(df, centroids)
        print(centroids)
        df = assignment(df, centroids)
        print(df)
        count += 1
        # Break if the clusters are the same after recalc the centroid
        if closest_centroids.equals(df['Cluster']):
            break
        # Break if it has gone through reassignment 20 times
        if count == 20:
            break
    # obtain final list of cluster
    final_list = []
    for i in np.unique(df.Cluster):
        print(i)
        print(np.array(df[df['Cluster'] == i]['Active']))
        final_list += [(np.array(df[df['Cluster'] == i]['Active']))]

    return final_list, df

def make_figure(df):
    """
    Will display the clusters

    Input: a dataframe with clusters
    Output: displays a 2D plot of the clusters, colored by the cluster column
    """

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    ax.scatter(df['principal component 1']
               , df['principal component 2']
               , c = df['Cluster']
               , s = 50)
    ax.grid()
    plt.show()
    return

# Add additional columns for each centroid calculating the distance from every point to that centroid
def assignment(df, centroids):
    """
    Will assign each active site to the nearest centroid

    Input: a dataframe and current centroids
    Output: new dataframe with updated cluster for each activestie
    """
    for i in range(len(centroids)):
        dist = np.array([])
        for j in range(len(df['Scores'])):
            temp = Euclidean_distance(df['Scores'][j], centroids[i])
            dist = np.append(dist, temp)
        df['distance_from_{}'.format(i)] = dist

    centroid_distance_cols = ['distance_from_{}'.format(i) for i in range(len(centroids))]
    # Selects the column of the centroid that they are closest too
    df['Cluster'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
    # remove distance from
    df['Cluster'] = df['Cluster'].map(lambda x: int(x.lstrip('distance_from_')))
    return df

# Get the mean of each newly generates cluster
def update(df, centroids):
    """
    # Get the mean of each newly generates cluster

    Input: a dataframe and current centroids
    Output: new centroids at the mean of each cluster
    """
    for i in range(len(centroids)):
        centroids[i] = np.average(np.array(df[df['Cluster'] == i]['Scores']).tolist(), axis=0)
    return centroids

def cluster_hierarchically(active_sites, k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances and K
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects) as well as a data frame describing the clusters
    """
    # Obtain similarity scores
    score_matrix = [[0]*len(active_sites) for i in  range(len(active_sites))]
    active = np.array([])
    X = -1
    for i in active_sites:
        X += 1
        Y = -1
        active = np.append(active, i.name)
        for j in active_sites:
            Y += 1
            score_matrix[X][Y] = compute_similarity(i,j)

    # Make score and cluster matrix
    d = {'Scores': score_matrix, 'Cluster': range(len(score_matrix)), 'Active':active}
    df = pd.DataFrame(data=d)
    print(df)
    # Use PCA to visulize all the vectors
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(score_matrix)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])
    df['principal component 1']=principalDf['principal component 1']
    df['principal component 2']=principalDf['principal component 2']
    dist = pd.DataFrame()
    for i in df.Cluster:
        d = np.array([])
        for j in df.Cluster:
            # gets the mean value of each cluster of points
            d = np.append(d, Euclidean_distance(np.average(np.array(df[df['Cluster'] == i]['Scores']).tolist(), axis=0),np.average(np.array(df[df['Cluster'] == j]['Scores']).tolist(), axis=0)))
        dist[i] = d
    # While cluster is less than K go through all the active sites and find the shortest active site distance and merge them together
    while k < len(np.unique(np.array([df.Cluster]))):
        print(dist)
        minval = float("inf")
        for i in np.unique(df.Cluster):
            col = (np.array(dist[i]))
            for j in np.unique(df.Cluster):
                if col[j] != 0 and minval > col[j]:
                        x = i
                        y = j
                        minval = col[j]
                        print(minval)
                        print(x)
                        print(y)
        # Cluster by having the mean dot
        mean = np.mean(np.array([(dist[x]),dist[y]]), axis =0)
        mean[x]=0
        print(mean)
        # Cluster by the max distance between each cluster
        MAX = np.fmax(np.array(dist[x]),np.array(dist[y]))
        MAX[x]=0
        MAX[y]=0
        # Clustered according to the max distance
        dist[x]=MAX
        df.Cluster=(df.replace(y,x)['Cluster']) #Replace the cluster of column y with x
    # obtain final list of cluster
    final_list = []
    for i in np.unique(df.Cluster):
        final_list += [(np.array(df[df['Cluster'] == i]['Active']))]

    return final_list, df


def Silhouette(df):
    # Has a problem if there is a cluster with just one


    # for k in kList:
    # df = cluster_by_partitioning(active_sites, 2)[1]
    score = []
    for c in np.unique(df.Cluster): # C is the current cluster being looked at
         same = (np.array(df[df['Cluster'] == c]['Scores'])) # Get all scores for that cluster
         # print(same)
         different = (np.array(df[df['Cluster'] != c]['Scores'])) # Get all scores for different cluster
         # print(different)
         for x in same: # For each point in the cluster of interest
            # print(x)
            same_list =[]
            dif_list = []
            for y in same: # Compare to each point in the same cluster excpet it's self
                # print(y)
                if not np.array_equal(x,y): # Makes sure not to compare to itself
                    same_list += [Euclidean_distance(x,y)] # Get distance of all dots in cluster in same list
            for y in different:
                dif_list += [Euclidean_distance(x,y)] # Get distance of all dots not in cluster in same list
            # print(same_list)
            # print(dif_list)
            mean = np.mean(same_list)/np.mean(dif_list) # Get the mean a(i)/b(i)
            score += [mean/np.max([np.mean(same_list), np.mean(dif_list)])] # Then divide the by the max to get the silouete score for that actibe site
    print('the final score is..')
    print(score)
    final = np.mean(score) # Finally calc the mean of all silouete scores
    print(final)

    return final

def find_optimal_score_P(active_sites):
    """
    Finds optimal clustering score using silhouette method for k-means clustering                                                            #

    Input: a list of ActiveSite instances
    Output: optimal K
    """
    m=0 #max
    k_opt=0 #optimal clustering
    for k in range(2,4):
        S = Silhouette(cluster_by_partitioning(active_sites, k)[1])
        print(S)
        if m < S:
            m=S
            k_opt=k
    print(k_opt)
    return k_opt

def find_optimal_score_H(active_sites):
    """
    Finds optimal clustering score using silhouette method for hierarchical clustering                                                            #

    Input: a list of ActiveSite instances
    Output: optimal K
    """
    m=0 #max
    k_opt=0 #optimal clustering
    for k in range(2,4):
        S = Silhouette(cluster_hierarchically(active_sites, k)[1])
        print(S)
        if m < S:
            m=S
            k_opt=k
    print(k_opt)
    return k_opt

def jaccard_index(df1,df2):
    # Compare how well two clustering algs performed
    cluster1 = np.array(df1.Cluster)
    cluster2 = np.array(df2.Cluster)
    total = len(df1.Cluster)
    same = 0
    for x in range(total):
        if cluster1[x] == cluster2[x]:
            same += 1

    print('the jaccard index is...')
    print(same/total)
    return same/total

    # Compare how well two clustering algs performed by taking the intersect over the union
    
    # psuedo code
    # The way I currently wrote it won't work because it doesn't account for clusters being the same, but named differently. Thus, I would need to find which clusters have two active sites in the same cluster compared to in different clusters as the scoring instead. To do this I could use a for loop comparing all the multiples of active sites to see that whether they are in the same list or not for each cluster. So that it is eventually total overlap/total active sites.    



