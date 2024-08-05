import numpy as np
import pandas as pd
import sys
from pyevtk.hl import pointsToVTK
from tqdm import tqdm
import argparse
import glob
from collections import defaultdict
import os
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
# sys.path.append('/work2/08264/baagee/frontera/cbgeopy/utils/')

def process_file(file_path):
    return pd.read_hdf(file_path, 'table')


def calculate_von_mises_stress(df):
    df['stress_magnitude'] = np.sqrt(
        0.5 * ((df['stress_xx'] - df['stress_yy']) ** 2 +
               (df['stress_yy'] - df['stress_zz']) ** 2 +
               (df['stress_zz'] - df['stress_xx']) ** 2) +
        3 * (df['tau_xy'] ** 2 + df['tau_yz'] ** 2 + df['tau_xz'] ** 2))


def calculate_displacement(df):
    df['displacement'] = np.sqrt(
        df['displacement_x'] ** 2 +
        df['displacement_y'] ** 2 +
        df['displacement_z'] ** 2)

def calculate_deviatoric_strain(df):
    df['deviatoric_strain'] = np.sqrt(
        (2 / 3) * ((df['strain_xx'] - df['strain_yy']) ** 2 +
                   (df['strain_yy'] - df['strain_zz']) ** 2 +
                   (df['strain_zz'] - df['strain_xx']) ** 2 +
                   3 * (df['gamma_xy'] ** 2 + df['gamma_yz'] ** 2 + df['gamma_xz'] ** 2)))


def write_vtk(df, save_path, vtk_type, t):
    x = np.ascontiguousarray(df['coord_x'])
    y = np.ascontiguousarray(df['coord_y'])
    z = np.ascontiguousarray(df['coord_z'])
    if x.size == 0:
        raise ValueError(f"`particles.h5` is empty at t={t}")

    if "stresses" in vtk_type:
        stress_xx = np.ascontiguousarray(df['stress_xx'])
        stress_yy = np.ascontiguousarray(df['stress_yy'])
        stress_zz = np.ascontiguousarray(df['stress_zz'])
        tau_xy = np.ascontiguousarray(df['tau_xy'])
        tau_yz = np.ascontiguousarray(df['tau_yz'])
        tau_xz = np.ascontiguousarray(df['tau_xz'])
        stress_magnitude = np.ascontiguousarray(df['stress_magnitude'])

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

    if "displacements" in vtk_type:
        displacement_x = np.ascontiguousarray(df['displacement_x'])
        displacement_y = np.ascontiguousarray(df['displacement_y'])
        displacement_z = np.ascontiguousarray(df['displacement_z'])
        displacement = np.ascontiguousarray(df['displacement'])

        pointsToVTK(
            f"{save_path}/particle_displacements-{t}",
            x, y, z,
            data={"displacement_x": displacement_x,
                  "displacement_y": displacement_y,
                  "displacement_z": displacement_z,
                  "displacement": displacement})

    if "strains" in vtk_type:
        strain_xx = np.ascontiguousarray(df['strain_xx'])
        strain_yy = np.ascontiguousarray(df['strain_yy'])
        strain_zz = np.ascontiguousarray(df['strain_zz'])
        gamma_xy = np.ascontiguousarray(df['gamma_xy'])
        gamma_yz = np.ascontiguousarray(df['gamma_yz'])
        gamma_xz = np.ascontiguousarray(df['gamma_xz'])
        deviatoric_strain = np.ascontiguousarray(df['deviatoric_strain'])

        pointsToVTK(
            f"{save_path}/particle_strains-{t}",
            x, y, z,
            data={"strain_xx": strain_xx,
                  "strain_yy": strain_yy,
                  "strain_zz": strain_zz,
                  "gamma_xy": gamma_xy,
                  "gamma_yz": gamma_yz,
                  "gamma_xz": gamma_xz,
                  "deviatoric_strain": deviatoric_strain})

    else:
        raise ValueError(f"VTK type option {vtk_type} does not exist")


def plot_stress_in_region_over_time(dfs, time_steps, x_range, y_range, z_range, stress_type='stress_magnitude',
                                    save_path=None):
    """
    Plot the stress over time for particles within a specified coordinate range.

    :param dfs: List of DataFrames, each containing particle data for a specific time step
    :param time_steps: List of time step values corresponding to each DataFrame
    :param x_range: Tuple of (min_x, max_x)
    :param y_range: Tuple of (min_y, max_y)
    :param z_range: Tuple of (min_z, max_z)
    :param stress_type: Type of stress to plot (default is 'stress_magnitude')
    :param save_path: Path to save the plot (if None, the plot will be displayed)
    """

    # Function to filter DataFrame based on coordinate range
    def filter_df(df):
        return df[
            (df['coord_x'] >= x_range[0]) & (df['coord_x'] <= x_range[1]) &
            (df['coord_y'] >= y_range[0]) & (df['coord_y'] <= y_range[1]) &
            (df['coord_z'] >= z_range[0]) & (df['coord_z'] <= z_range[1])
            ]

    # Calculate statistics for each time step
    stats = []
    for df in dfs:
        filtered_df = filter_df(df)
        stats.append({
            'mean': filtered_df[stress_type].mean(),
            'median': filtered_df[stress_type].median(),
            'min': filtered_df[stress_type].min(),
            'max': filtered_df[stress_type].max(),
            'std': filtered_df[stress_type].std()
        })

    stats_df = pd.DataFrame(stats, index=time_steps)

    # Create the plot
    plt.figure(figsize=(12, 8))

    plt.plot(stats_df.index, stats_df['mean'], label='Mean', color='blue')
    plt.fill_between(stats_df.index, stats_df['mean'] - stats_df['std'],
                     stats_df['mean'] + stats_df['std'], alpha=0.2, color='blue')

    plt.plot(stats_df.index, stats_df['median'], label='Median', color='green')
    plt.plot(stats_df.index, stats_df['min'], label='Min', color='red')
    plt.plot(stats_df.index, stats_df['max'], label='Max', color='purple')

    plt.xlabel('Time Step')
    plt.ylabel(stress_type.replace('_', ' ').title())
    plt.title(f'{stress_type.replace("_", " ").title()} Over Time in Specified Region')
    plt.legend()
    plt.grid(True)

    # Add text box with coordinate ranges
    text_box = f'X range: {x_range}\nY range: {y_range}\nZ range: {z_range}'
    plt.text(0.95, 0.05, text_box, transform=plt.gca().transAxes, fontsize=9,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    # Save or show the plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    else:
        plt.show()



def process_time_step(t, h5s_at_t, result_dir, save_path, vtk_type):
    dfs = []
    for h5_name in h5s_at_t:
        h5_path = os.path.join(result_dir, h5_name)
        df = process_file(h5_path)
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    if "stresses" in vtk_type:
        calculate_von_mises_stress(df)
    if "displacements" in vtk_type:
        calculate_displacement(df)
    if "strains" in vtk_type:
        calculate_deviatoric_strain(df)

    write_vtk(df, save_path, vtk_type, t)

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hdf5 trajectories to vtk.')
    parser.add_argument('--result_dir',
                        default='/scratch1/08264/baagee/cbgeopy-scratch/simulations/fundao3d-8-2/results/sand3d-1/',
                        help="Path to hdf5 files to consume.")
    parser.add_argument('--vtk', type=str,
                        default=['stresses', 'strains'],
                        nargs='+',
                        help="Options to vtk return: ['displacements', 'stresses', 'strains']")
    parser.add_argument('--max_time',
                        default=200000,
                        help="Max timestep to convert")
    parser.add_argument('--mpi', type=int,
                        default=32,
                        help="Number of MPI used for mpm runtime")
    parser.add_argument('--write_workers', type=int,
                        default=32,
                        help="Number of workers for vtk writing")
    parser.add_argument('--save_path',
                        default='/scratch1/08264/baagee/cbgeopy-scratch/simulations/fundao3d-8-2/results/sand3d-1/',
                        help="Save path")
    args = parser.parse_args()
    stress_plot = True

    result_dir = args.result_dir
    max_time = args.max_time
    vtk = args.vtk
    mpi = args.mpi
    save_path = args.save_path
    write_workers = args.write_workers

    file_paths = glob.glob(f'{result_dir}/particles*.h5')
    files_names = [file_path.split("/")[-1] for file_path in file_paths]

    # Initialize a dictionary to hold the groups
    files_by_time = defaultdict(list)

    # Group the file names by the same timestep
    for files_name in files_names:
        # Extract the time component from the filename
        time = files_name.split('-')[-1].split('.')[0]
        if max_time is not None:
            if int(time) <= int(max_time):
                # Add the file to the corresponding time group
                files_by_time[time].append(files_name)
        else:
            files_by_time[time].append(files_name)

    # Reorder the dictionary by timestep ascending order
    files_by_time = dict(sorted(files_by_time.items(), key=lambda item: int(item[0])))

    all_dfs = []
    time_steps = []

    with ProcessPoolExecutor(max_workers=write_workers) as executor:
        futures = [
            executor.submit(process_time_step, t, h5s_at_t, result_dir, save_path, vtk)
            for t, h5s_at_t in files_by_time.items()
        ]

        for future, t in zip(tqdm(futures), files_by_time.keys()):
            df = future.result()
            all_dfs.append(df)
            time_steps.append(int(t))

    # Example usage of the new plotting function
    if stress_plot:
        plot_stress_in_region_over_time(
            all_dfs,
            time_steps,
            x_range=(800, 810),  # Replace with your desired x range
            y_range=(400, 410),  # Replace with your desired y range
            z_range=(740, 750),  # Replace with your desired z range
            stress_type='stress_zz',
            save_path=f'{save_path}/stress_in_region_over_time.png'
        )
