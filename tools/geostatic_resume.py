import numpy as np
import pandas as pd
import os
import json
# import subprocess


# Number of MPI tasks
mpi = 32
result_dir = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/results/'
uuid_stress_equilibrium = 'sand3d-le'
uuid_save = 'sand3d-resume'
timestep_undeform = '0000000'
timestep_stress_equilibrium = '0200000'
results_stress_equilibrium = f'{result_dir}/{uuid_stress_equilibrium}/'
entity_set_path = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/entity_sets.json'


# Get entity set file
with open(entity_set_path) as f:
    particle_sets = json.load(f)['particle_sets']


def get_h5(uuid_dir, timestep):
    # Create an empty list to store DataFrames
    dfs = []

    # Iterate over different files from MPI and append to list
    for i in range(mpi):
        file = f'particles-{i}_{mpi}-{timestep}.h5'  # ex) particles-26_32-0120000.h5
        h5_path = os.path.join(uuid_dir, file)
        dfs.append(pd.read_hdf(h5_path, 'table'))

    # Concatenate all DataFrames
    df = pd.concat(dfs, ignore_index=True)

    return df

# # df by paritcles
# dfp0 = df.loc[df['material_id'] == 1]
# dfp1 = df.loc[df['material_id'] == 0]

df_undeformed = get_h5(
    f'{result_dir}/{uuid_stress_equilibrium}', timestep_undeform)
df_stress_equilibrium = get_h5(
    f'{result_dir}/{uuid_stress_equilibrium}', timestep_stress_equilibrium)

# Copy the df when no deformation.
df_geostatic = df_undeformed.copy()

# Define stress columns
stress_columns = list(df_geostatic.columns[16:28])

# Merge dataframes on 'id'
merged_df = df_geostatic[['id']].merge(
    df_stress_equilibrium[['id'] + stress_columns], on='id', how='left')

# Update the stress values in df_b
df_geostatic[stress_columns] = merged_df[stress_columns]

if not os.path.exists(f"{result_dir}/{uuid_save}/"):
    os.makedirs(f"{result_dir}/{uuid_save}/", exist_ok=True)
df_geostatic.to_csv(f"{result_dir}/{uuid_save}/particles{timestep_undeform}.csv", index=False)




# # Insert the deformed stress to the stabilized stress.
# for particle_set in particle_sets:
#     id = particle_set['id']
#     pset = particle_set['set']
#     pset = [0, 1, 2]
#
#     # Copy the stresses at the stress equilibrium to the none-deformed original particles.
#     mask = df_stress_equilibrium['id'].isin(pset)
#     stress_values = df_stress_equilibrium.loc[mask, stress_columns]
#     df_geostatic.loc[mask, stress_columns] = stress_values.values


a=1