import os
import pandas as pd


mpi = 32
timestep = '0070000'
result_dir = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/fundao3d-8-2/results/sand3d-1/'

# Create an empty list to store DataFrames
dfs = []

# Iterate over different files from MPI and append to list
for i in range(mpi):
    file = f'particles-{i}_{mpi}-{timestep}.h5'  # ex) particles-26_32-0120000.h5
    h5_path = os.path.join(result_dir, file)
    dfs.append(pd.read_hdf(h5_path, 'table'))

# Concatenate all DataFrames
df = pd.concat(dfs, ignore_index=True)
a=1