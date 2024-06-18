import numpy as np
import pandas as pd
from pyevtk.hl import pointsToVTK
from tqdm import tqdm
from utils import get_h5
import argparse
import glob
from collections import defaultdict
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hdf5 trajectories to vtk.')
    parser.add_argument('--result_dir', default='/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/results/sand3d-resume/', help="Path to hdf5 files to consume.")
    parser.add_argument('--vtk', default='stresses', help="Options to vtk return: ['displacements', 'stresses', 'strains']")
    parser.add_argument('--mpi', default='32', help="Number of MPI")
    parser.add_argument('--save_path', default='/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand_layers_geostatic/results/', help="Save path")
    args = parser.parse_args()

    result_dir = args.result_dir
    vtk = args.vtk
    mpi = args.mpi
    save_path = args.save_path

    file_paths = glob.glob(f'{result_dir}/particles*.h5')
    files_names = [file_path.split("/")[-1] for file_path in file_paths]

    # Initialize a dictionary to hold the groups
    files_by_time = defaultdict(list)

    # Group the file names by the same timestep
    for files_name in files_names:
        # Extract the time component from the filename
        time = files_name.split('-')[-1].split('.')[0]
        # Add the file to the corresponding time group
        files_by_time[time].append(files_name)

    # Reorder the dictionary by timestep ascending order
    files_by_time = dict(sorted(files_by_time.items(), key=lambda item: int(item[0])))

    for t, h5s_at_t in tqdm(files_by_time.items()):
        # df at time t
        dfs = []
        for h5_name in h5s_at_t:
            h5_path = os.path.join(result_dir, h5_name)
            dfs.append(pd.read_hdf(h5_path, 'table'))

        # Concatenate all DataFrames
        df = pd.concat(dfs, ignore_index=True)

        # Extract positions
        x = np.ascontiguousarray(df['coord_x'])
        y = np.ascontiguousarray(df['coord_y'])
        z = np.ascontiguousarray(df['coord_z'])

        # Save vtp
        if vtk == "stresses":
            # Calculate von Mises stress magnitude
            df['stress_magnitude'] = np.sqrt(
                0.5 * ((df['stress_xx'] - df['stress_yy']) ** 2 +
                       (df['stress_yy'] - df['stress_zz']) ** 2 +
                       (df['stress_zz'] - df['stress_xx']) ** 2) +
                3 * (df['tau_xy'] ** 2 + df['tau_yz'] ** 2 + df['tau_xz'] ** 2))


            # Extract stress components and magnitude
            stress_xx = np.ascontiguousarray(df['stress_xx'])
            stress_yy = np.ascontiguousarray(df['stress_yy'])
            stress_zz = np.ascontiguousarray(df['stress_zz'])
            tau_xy = np.ascontiguousarray(df['tau_xy'])
            tau_yz = np.ascontiguousarray(df['tau_yz'])
            tau_xz = np.ascontiguousarray(df['tau_xz'])
            stress_magnitude = np.ascontiguousarray(df['stress_magnitude'])

            # Write to VTP
            pointsToVTK(
                f"{save_path}/particle_stresses-{t}",
                x, y, z,
                data={"stress_xx": stress_xx,
                      "stress_yy": stress_yy,
                      "stress_zz": stress_zz,
                      "tau_xy": tau_xy,
                      "tau_yz": tau_yz,
                      "tau_xz": tau_xz,
                      "stress_magnitude": stress_magnitude})

        elif vtk == "displacements":
            # Calculate total displacement
            df['displacement'] = np.sqrt(
                df['displacement_x'] ** 2 +
                df['displacement_y'] ** 2 +
                df['displacement_z'] ** 2)

            # Extract displacement components
            displacement_x = np.ascontiguousarray(df['displacement_x'])
            displacement_y = np.ascontiguousarray(df['displacement_y'])
            displacement_z = np.ascontiguousarray(df['displacement_z'])
            displacement = np.ascontiguousarray(df['displacement'])

            # Write to VTP
            pointsToVTK(
                f"{save_path}/particle_displacements-{t}",
                x, y, z,
                data={"displacement_x": displacement_x,
                      "displacement_y": displacement_y,
                      "displacement_z": displacement_z,
                      "displacement": displacement})

