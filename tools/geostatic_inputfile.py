### Note ###
# This code gathers distributed `.h5` files at the estimated stress equilibrium timestep, and
# creates `particles_stresses.txt` so that we can start simulation with this stress state.
############

import numpy as np
import pandas as pd
import os
import json
import argparse


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


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Gather .h5 files and generate particles_stresses.txt")

    parser.add_argument(
        '--mpi', type=int,
        default=4,
        help="Number of MPI processes")
    parser.add_argument(
        '--result_dir', type=str,
        default="./",
        help="Directory for simulation results")
    parser.add_argument(
        '--result_subdir', type=str,
        default="results/sand2d/",
        help="UUID for stress equilibrium")
    parser.add_argument(
        '--timestep_undeform', type=str,
        default="000000",
        help="Timestep for undeformed state")
    parser.add_argument(
        '--timestep_stress_equilibrium', type=str,
        default=None,
        help="Timestep for stress equilibrium")
    parser.add_argument(
        '--k0', type=float,
        default=None,
        help="Compute stress based on the height of the mass (i.e., k0 * rho * g * h)."
             "This only works for rectangular mass")
    parser.add_argument(
        '--unit_weight', type=float,
        default=None,
        help="Compute stress based on the height of the mass (i.e., k0 * rho * g * h)."
             "This only works for rectangular mass")

    args = parser.parse_args()

    # Assign inputs from arguments
    mpi = args.mpi
    result_dir = args.result_dir
    result_subdir = args.result_subdir
    timestep_undeform = args.timestep_undeform
    timestep_stress_equilibrium = args.timestep_stress_equilibrium
    k0 = args.k0
    unit_weight = args.unit_weight

    # # Debug
    # mpi = 4
    # result_dir = "/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_multilayers/sim-1/"
    # result_subdir = "results/sand2d-le/"
    # timestep_undeform = "0000"
    # # timestep_stress_equilibrium = "02250"
    # k0 = 0.5
    # unit_weight = 1800

    df_undeformed = get_h5(
        f'{result_dir}/{result_subdir}', timestep_undeform, n_mpis=mpi)

    if all([k0, unit_weight]):
        if timestep_stress_equilibrium:
            print("Waring: `timestep_stress_equilibrium` is not necessary in k0 mode")

        # get particle positions
        particles = df_undeformed[['coord_x', 'coord_y', 'coord_z']].to_numpy()
        # compute k0 stress based
        z_max = particles[:, 1].max()
        particle_stress = np.zeros(
            (np.shape(particles)[0], 6))  # second axis is for sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_zx,
        particle_stress[:, 0] = k0 * (z_max - particles[:, 1]) * unit_weight  # K0*H*Unit_Weight
        particle_stress[:, 1] = (z_max - particles[:, 1]) * unit_weight  # H*Unit_Weight
        particle_stress[:, 2] = np.zeros(np.shape(particles)[0])
        particle_stress[:, 3] = np.zeros(np.shape(particles)[0])
        particle_stress[:, 4] = np.zeros(np.shape(particles)[0])
        particle_stress[:, 5] = np.zeros(np.shape(particles)[0])
        # Create a header string with the number of particles
        header = str(len(particle_stress))

        # Save the array to a text file
        np.savetxt(
            f"{result_dir}/particles-stresses.txt",
            particle_stress,
            delimiter="\t", fmt="%.4f", header=header, comments="")

    else:
        if timestep_stress_equilibrium is None:
            raise ValueError("`timestep_stress_equilibrium` should be passed")

        df_stress_equilibrium = get_h5(
            f'{result_dir}/{result_subdir}', timestep_stress_equilibrium, n_mpis=mpi
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


if __name__ == "__main__":
    main()