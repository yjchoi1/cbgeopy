import os
import pandas as pd


mpi = 4
timestep = '0080000'
result_dir = '/work2/08264/baagee/frontera/cbgeopy/simulations/fundao-2d/results/sand2d/'

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