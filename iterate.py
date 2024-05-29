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



