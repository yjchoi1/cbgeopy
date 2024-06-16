import numpy as np
import pandas as pd
import os
import json



# Number of MPI tasks
mpi = 32
result_dir = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/results/'
uuid_stress_equilibrium = 'sand3d-le'
uuid_save = 'sand3d-resume'
timestep_undeform = '0000000'
timestep_stress_equilibrium = '0200000'
results_stress_equilibrium = f'{result_dir}/{uuid_stress_equilibrium}/'
entity_set_path = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/entity_sets.json'


def get_h5(uuid_dir, timestep, n_mpis):
    # Create an empty list to store DataFrames
    dfs = []

    # Iterate over different files from MPI and append to list
    for i in range(n_mpis):
        file = f'particles-{i}_{n_mpis}-{timestep}.h5'  # ex) particles-26_32-0120000.h5
        h5_path = os.path.join(uuid_dir, file)
        dfs.append(pd.read_hdf(h5_path, 'table'))

        # Concatenate all DataFrames
        df = pd.concat(dfs, ignore_index=True)

    return df

df_undeformed = get_h5(
    f'{result_dir}/{uuid_stress_equilibrium}', timestep_undeform, n_mpis=mpi)
df_stress_equilibrium = get_h5(
    f'{result_dir}/{uuid_stress_equilibrium}', timestep_stress_equilibrium, n_mpis=mpi
)
# Copy the df when no deformation.
df_geostatic = df_undeformed.copy()


# Define stress columns
stress_columns = list(df_geostatic.columns[16:16+6])

# Merge `df_stress_equilibrium` based on 'id' of `df_geostatic` which follows `df_undeformed`
merged_df = df_geostatic[['id']].merge(
    df_stress_equilibrium[['id'] + stress_columns], on='id', how='left')

# Update the stress values in `df_geostatic`
df_geostatic[stress_columns] = merged_df[stress_columns]
df_geostatic_sorted = df_geostatic.sort_values(by='id').reset_index(drop=True)

# Write the number of stressed particles
with open(f"{result_dir}/particles-stresses.txt", "w") as f:
    f.write(f"{len(df_stress_equilibrium)} \n")
# Append the dataframe to the file with tab separation
df_geostatic_sorted[stress_columns].to_csv(
    f"{result_dir}/particles-stresses.txt", mode='a', sep='\t', index=False,  header=False)
