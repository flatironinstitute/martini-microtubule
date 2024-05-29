import numpy as np
import networkx as nx
import pickle
import os
from shutil import copy2, copytree, rmtree
from MDAnalysis import Universe

# Constants and Configuration
ALPHA = 1050
RC = 0.9
COEFF = (ALPHA / (RC**4)) * 2.5206
CURRENT_ITER = 2
PREV_ITER = CURRENT_ITER - 1
NEXT_ITER = CURRENT_ITER + 1
next_iter_path = f'../iter{NEXT_ITER}'

# File Operations
def safe_copy(src, dst, directory=False):
    """Copy files or directories safely."""
    if directory:
        if os.path.exists(dst):
            rmtree(dst)
        copytree(src, dst)
    else:
        copy2(src, dst)

# Securely copying directories and files
safe_copy(f'../iter{PREV_ITER}/Bag/', './Bag/', directory=True)
safe_copy(f'../iter{PREV_ITER}/attributes/', './attributes/', directory=True)
safe_copy('./attributes/enm', './enm_cg_prev')
safe_copy('./attributes/enm_a1b2', './enm_a1b2_cg_prev')

# Function to concatenate files safely
def concatenate_files(output_file, *input_files):
    with open(output_file, 'wb') as outfile:
        for fname in input_files:
            with open(fname, 'rb') as infile:
                outfile.write(infile.read())

concatenate_files('./cg_prev', './enm_cg_prev', './enm_a1b2_cg_prev')

# Loading and configuring MDAnalysis
univ = Universe('protofilament.gro')
atom_selections = {
    'b1': (0, 1022, 424),
    'a1': (1022, 2067, 432),
    'b2': (2067, 2067+1022, 424),
    'a2': (2067+1022, 2067+1022+1045, 432)
}
atoms = {key: univ.atoms[start:end].select_atoms('name BB')[:count] for key, (start, end, count) in atom_selections.items()}
TEMP = sum(atoms.values(), start=next(iter(atoms.values())))

# Loading previous iteration for edge simulation
u = Universe(f'../iter{PREV_ITER}/cg_iter1_pbc.tpr', f'../iter{PREV_ITER}/cg_iter1_pbc.xtc')
PROTCG = [u.atoms[2067*i:2067*(i+1)] for i in range(12)]
prev_atoms = {key: PROTCG[idx][:1022].select_atoms('name BB')[:424] if 'beta' in key else PROTCG[idx][1022:].select_atoms('name BB')[:432] 
              for key, idx in {'beta1': 5, 'alpha1': 5, 'beta2': 6, 'alpha2': 6}.items()}
PREV = sum(prev_atoms.values(), start=next(iter(prev_atoms.values())))

# Read pickle data and compose graphs
def read_pickle(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

graph_files = ['./attributes/b1b1', './attributes/a1a1', './attributes/a1b1', './attributes/a1b2']
graphs4 = [read_pickle(f) for f in graph_files]
graphs3 = [read_pickle(f) for f in graph_files[:3]]
G_set4 = nx.compose_all(graphs4)
G_set3 = nx.compose_all(graphs3)

# Loading previous CG data
CG_prev_enm = np.loadtxt('./cg_prev')

#From all-atom simulations ...

##AA_ID_MATCH; Creates G_aa from atomisic simulation
G_aa = nx.Graph()
List_of_Edges = list(G_set4.edges)
for ed in List_of_Edges:
    ei = int(ed[0][:-3])
    ej = int(ed[1][:-3])
    dmean = G_set4.edges[ed[0], ed[1]]['weight'][0]/10.0
    variance = G_set4.edges[ed[0], ed[1]]['weight'][1]
    if ed[0][-2:] == 'a1':
        atm1 = a1[ei]
    if ed[0][-2:] == 'b1':
        atm1 = b1[ei]
    if ed[0][-2:] == 'b2':
        atm1 = b2[ei]
    if ed[1][-2:] == 'a1':
        atm2 = a1[ej]
    if ed[1][-2:] == 'b1':
        atm2 = b1[ej]
    if ed[1][-2:] == 'b2':
        atm2 = b2[ej]
    G_aa.add_edge(str(atm1.id), str(atm2.id), weight=[dmean, variance])
    
    

def calculate_trajectory_variances(u, atom_pairs, num_frames):
    """
    Calculate variances of distances between pairs of atoms over specified trajectory frames.
    
    Args:
        u (Universe): MDAnalysis Universe object containing trajectory data.
        atom_pairs (list): List of tuples, each containing two Atom objects to measure distances.
        num_frames (int): Number of last frames to analyze.
    
    Returns:
        np.ndarray: Variances of the distances between atom pairs across the given frames.
    """
    distance_data = []
    for ts in u.trajectory[-num_frames:]:  # Consider only the last num_frames
        distances = [np.linalg.norm(atom1.positions - atom2.positions) for atom1, atom2 in atom_pairs]
        distance_data.append(distances)
    return np.var(np.array(distance_data), axis=0)

def construct_previous_cg_graph(CG_data, variances):
    """
    Construct a graph from previous CG data with variance as weights.
    
    Args:
        CG_data (np.ndarray): Array of previous CG data with edge indices and other properties.
        variances (np.ndarray): Array of variance data for the edges.
    
    Returns:
        nx.Graph: A NetworkX graph with edges and associated variance weights.
    """
    previous_cg = nx.Graph()
    for (index_data, variance) in zip(CG_data, variances):
        # index_data expected format: [index_i, index_j, ..., property_k]
        node1, node2, weight = int(index_data[0]), int(index_data[1]), index_data[4]
        previous_cg.add_edge(str(node1), str(node2), weight=[weight, variance])
    return previous_cg

# Setup variables
offset = beta1[0].id - 1  # Adjust offset based on atom ID
ECG_indices_i = CG_prev_enm[:, 0].astype(int) + offset
ECG_indices_j = CG_prev_enm[:, 1].astype(int) + offset
atom_pairs = [(PREV.select_atoms(f'id {idx_i}')[0], PREV.select_atoms(f'id {idx_j}')[0])
              for idx_i, idx_j in zip(ECG_indices_i, ECG_indices_j)]

# Calculate variances and create graph
variances = calculate_trajectory_variances(u, atom_pairs, 1000)
Previous_CG = construct_previous_cg_graph(CG_prev_enm, variances)


def update_cg_model(G_aa, Previous_CG, coeff):
    """
    Updates the new CG model based on differences in variance and previous CG data.
    
    Args:
        G_aa (nx.Graph): Graph of atomistic simulation data.
        Previous_CG (nx.Graph): Graph of previous CG data.
        coeff (float): Coefficient used in the force constant calculation.
    
    Returns:
        nx.Graph: Updated new CG model.
    """
    GCG_new = nx.Graph()
    edge_list = np.array(list(G_aa.edges))
    ELi, ELj = edge_list[:, 0].astype(int), edge_list[:, 1].astype(int)

    with open('enm', 'w+') as enm, open('enm_a1b2', 'w+') as enm_a1b2:
        for ei, ej in zip(ELi, ELj):
            ai, aj = str(ei), str(ej)
            try:
                var_cg = Previous_CG.edges[ai, aj]['weight'][1] / 100.0
                var_aa = G_aa.edges[ai, aj]['weight'][1] / 100.0
                d = G_aa.edges[ai, aj]['weight'][0]
                k_prev = Previous_CG.edges[ai, aj]['weight'][0]
                Delta = var_aa - var_cg
                k_trial = k_prev - (coeff * Delta)
                k_new = max(0.0, k_trial)  # Ensure that k_new is not negative

                # Select the file based on edge indices
                file = enm_a1b2 if ei > 2067 or ej > 2067 else enm
                file.write(f"{ei}\t{ej}\t1\t{d}\t{k_new}\n")
                GCG_new.add_edge(ei, ej, weight=[d, k_new, k_prev])
            except KeyError as e:
                print(f"Key error: {e} - Possible missing edge in Previous_CG.")

    return GCG_new


# Assuming G_aa and Previous_CG are already defined and populated
GCG_new = update_cg_model(G_aa, Previous_CG, COEFF)

def create_directory(directory_path):
    """Safely create directory if it doesn't already exist."""
    os.makedirs(directory_path, exist_ok=True)



# Copying ENM files to attributes and Bag directories
safe_copy('enm*', './attributes/')
safe_copy('enm*', './Bag/')

# Create next iteration directory and copy 'iterate.py'
create_directory(next_iter_path)
shutil.copy('iterate.py', f'{next_iter_path}/')
