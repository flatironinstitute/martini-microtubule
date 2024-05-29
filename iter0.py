import numpy as np
from MDAnalysis import *
import os
from MDAnalysis.analysis.distances import distance_array
import networkx as nx
import pickle



def outliers_modified_z_score(ys, ref, threshold=3.5):
    """
    Determine if a reference value is an outlier based on the modified Z-score method.

    Parameters:
    - ys (list or array): Input values to compare against.
    - ref (float): Reference value to test as an outlier.
    - threshold (float): Threshold for the Z-score above which the value is considered an outlier.

    Returns:
    - bool: True if the reference value is not an outlier, False otherwise.
    """
    ys_arr = np.array(ys)
    median_y = np.median(ys_arr)
    deviations = np.abs(ys_arr - median_y)
    mad_y = np.median(deviations)  # Median absolute deviation

    if mad_y == 0:  # Avoid division by zero
        return False  # Assuming that all ys are the same, and thus ref is not an outlier

    modified_z_scores = 0.6745 * (ref - median_y) / mad_y
    return np.abs(modified_z_scores) < threshold




def create_graph(unit, set1=beta1, set2=beta1, ref='b1b1'):
    # Calculate reference distances for the first trajectory
    ref_distances = distance_array(set1.positions, set2.positions)
    
    # Calculate distances for the last 500 frames
    distances = [distance_array(set1.positions, set2.positions) for ts in unit.trajectory[-500:]]
    distances = np.array(distances)
    
    # Compute mean, variance, and standard deviation of distances
    mean_distances = np.mean(distances, axis=0)
    variance_distances = np.var(distances, axis=0)
    
    # Create a graph
    graph = nx.Graph()
    close_pairs = np.where((mean_distances < 9.0) & (mean_distances >= 5.0))
    
    # Add edges to the graph based on conditions
    for i, j in zip(*close_pairs):
        if ref[:2] == ref[2:] and np.abs(i - j) > 2:
            condition = outliers_modified_z_score(distances[:, i, j], ref_distances[i, j], threshold=3.5)
        else:
            condition = outliers_modified_z_score(distances[:, i, j], ref_distances[i, j], threshold=3.5)        
        if condition:
            graph.add_edge(f"{i}_{ref[:2]}", f"{j}_{ref[2:]}", weight=[mean_distances[i, j], variance_distances[i, j]])
    
    # Save graph to file
    with open(f'./attributes/{ref}', 'wb') as f:
        pickle.dump(graph, f)
    
    return graph


# Load the universe with topology and trajectory files
u = Universe('/mnt/home/asahoo/ceph/setups/MicroT/7sja/6DPV/Z3/6dpv_pbc.tpr',
             '/mnt/home/asahoo/ceph/setups/MicroT/7sja/6DPV/Z3/6dpv_pbc_2.xtc')

# Select protein atoms and divide into segments for processing
prot = u.select_atoms('protein')
PROTAA = [prot[13774*i:13774*(i+1)] for i in range(12)]

# Define subsets of atoms for different chains and types
beta1 = PROTAA[6][:6840].select_atoms('name CA')[:424]
alpha1 = PROTAA[6][6840:].select_atoms('name CA')[:432]
beta2 = PROTAA[7][:6840].select_atoms('name CA')[:424]
alpha2 = PROTAA[7][6840:].select_atoms('name CA')[:432]


#Check for chainID nums
#beta-beta:Y
#alpha-alpha: W
#intra: aY, WU
#inter: YW
#CG
# Load the universe from the protofilament file
univ = Universe('protofilament.gro')

# Select specific atoms groups by their names and indices
b1 = univ.atoms[:1022].select_atoms('name BB')[:424]
a1 = univ.atoms[1022:2067].select_atoms('name BB')[:432]
b2 = univ.atoms[2067:2067 + 1022].select_atoms('name BB')[:424]
a2 = univ.atoms[2067 + 1022:2067 + 1022 + 1045].select_atoms('name BB')[:432]

# Definitions for easy reference in creating graphs
ref_labels = {
    'b1b1': (b1, b1),
    'b2b2': (b2, b2),
    'a1a1': (a1, a1),
    'a2a2': (a2, a2),
    'a1b1': (a1, b1),
    'a2b2': (a2, b2),
    'a1b2': (a1, b2)
}


# Create graphs for each interaction pair
graphs = {label: Create_Graph(univ, set1=sets[0], set2=sets[1], ref=label) for label, sets in ref_labels.items()}

# Combine the graphs
G_set4 = nx.compose_all(graphs.values())

# Calculate the minimum variance from all edges in the final composed graph
Variances = [G_set4.edges[edge]['weight'][1] for edge in G_set4.edges]
Var_min = np.min(Variances)

# Generate ENM files for subsets of the composed graph and specific graphs
gen_enm(nx.compose_all([graphs['b1b1'], graphs['a1a1'], graphs['a1b1']]), Var_min, ENM='enm', ENM_COMPARE='enm_compare')
gen_enm(graphs['a1b2'], Var_min, ENM='enm_a1b2', ENM_COMPARE='enm_compare_a1b2')



'''
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/enm /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/enm_a1b2 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/enm_compare /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/enm_compare_a1b2 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')


os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/a1a1 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/b1b1 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/a1b1 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')
os.system('cp /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/attributes/a1b2 /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Bag/')

os.system('python /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/itp_modify_withallrestrained.py')
os.system('python /mnt/home/asahoo/ceph/setups/MicroT/7sja/mart_6DPV/makeCG1/Iterative/iter1/Create/itp_modify_withallrestrained_ends.py')
'''


