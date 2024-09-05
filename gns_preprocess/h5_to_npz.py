from tqdm import tqdm
import glob
import json
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple


def scale_positions(positions, origin, scaling_factor):
    """
    Scales the positions of particles with respect to a specified origin by a specified factor.

    Parameters:
    - positions: np.ndarray of shape (n_timesteps, n_particles, 2),
                 the positions of particles over time in 2D space.
    - origin: np.ndarray of shape (2,), the origin with respect to which the scaling is performed.
    - scaling_factor: float, the factor by which to scale the positions.

    Returns:
    - new_positions: np.ndarray of the same shape as positions, the scaled positions.
    """
    # Step 1: Subtract the origin from the positions
    shifted_positions = positions - origin

    # Step 2: Apply the scaling factor
    scaled_positions = shifted_positions * scaling_factor

    # Step 3: Add the origin back to the positions
    new_positions = scaled_positions + origin

    return new_positions


def concat_h5s_at_t(h5s_at_t: List[str], result_dir: str) -> pd.DataFrame:
    dfs = []
    columns = ['id', 'coord_x', 'coord_y', 'coord_z', 'material_id']
    for h5_name in h5s_at_t:
        h5_path = os.path.join(result_dir, h5_name)
        df = pd.read_hdf(h5_path, 'table', columns=columns)
        dfs.append(df)

    # For the mpm result from mpi job, the particle order might be different for each timestep.
    # so, use sort to reorder the data based on particle id with ascending order
    df_concat = pd.concat(dfs, ignore_index=True).sort_values(by='id', ascending=True)
    return df_concat


def process_time_step(t, h5s_at_t, result_dir, n_dims):
    df = concat_h5s_at_t(h5s_at_t, result_dir)

    # positions
    positions = df[['coord_x', 'coord_y', 'coord_z'][:n_dims]].to_numpy()

    return t, positions, df


def convert_hd5_to_npz(
        n_dims: int,
        result_dir: str,
        mpm_input_path: str,
        save_path: str,
        use_material_feature: bool,
        scale: Dict = None,
        max_time: int = None,
        mpi: int = 1,
        dt: float = 1.0):
    """
    Read the results of CB-Geo MPM and convert it to npz format, which is compatible with GNS.
    It can
    Args:
        n_dims ():
        result_dir ():
        mpm_input_path ():
        save_path ():
        use_material_feature ():
        scale (Dict): {"origin": [x, y, z], "factor": 2.0}
        max_time ():
        mpi ():
        dt ():

    Returns:

    """
    file_paths = glob.glob(f'{result_dir}/particles*.h5')
    files_names = [file_path.split("/")[-1] for file_path in file_paths]

    # Initialize a dictionary to hold the groups
    files_by_time = defaultdict(list)

    # Group the file names by the same timestep
    # TODO: it can do both mpi and serial format
    #   For example, `particles-000000.h5` or `particles-0_8-001000.h5
    #   Check if it can do.
    for files_name in files_names:
        # Extract the time component from the filename
        time = files_name.split('particles')[-1].split('-')[-1].split('.')[0]
        if max_time is not None:
            if int(time) <= int(max_time):
                # Add the file to the corresponding time group
                files_by_time[time].append(files_name)
        else:
            files_by_time[time].append(files_name)

    # Reorder the dictionary by timestep ascending order
    files_by_time = dict(sorted(files_by_time.items(), key=lambda item: int(item[0])))

    n_particles_reference = None
    positions_over_time = []

    # Parallelize the processing of each time step
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(
            process_time_step, t, h5s_at_t, result_dir, n_dims)
                   for t, h5s_at_t in files_by_time.items()]
        for future in tqdm(as_completed(futures)):
            t, positions, df = future.result()
            # Store positions over time
            positions_over_time.append((t, positions))

            # print(f"Processing timestep {t}")

            # Check if the particle disappears
            n_particles = len(positions)
            if n_particles_reference is None:
                n_particles_reference = n_particles
            else:
                if n_particles != n_particles_reference:
                    raise ValueError(f"Number of particles at timestep {t} is not the same as previous one.")

    positions_over_time.sort(key=lambda x: x[0])
    positions_list = [positions for _, positions in positions_over_time]

    # concat positions over time
    positions_over_time = np.stack(
        positions_list, axis=0)[::scale["timestep_stride"]]

    positions_over_time = scale_positions(
        positions_over_time, scale["origin"], scale["factor"])

    fig, ax = plt.subplots()
    plottime = -1
    ax.scatter(positions_over_time[plottime, :, 0], positions_over_time[plottime, :, 1])
    plt.show()

    # material feature
    if use_material_feature:
        # open mpm.json
        with open(mpm_input_path, 'r') as file:
            mpm_json = json.load(file)

        # get material id and associated properties
        materials = mpm_json["materials"]

        # get the map between material_id and corresponding property
        particle_material_mapping = {}  # {material_id: phi_deg}
        for material in materials:
            particle_material_mapping[material["id"]] = material["friction"]

        # make material_feature based on the id.
        material_feature = df['material_id'].map(particle_material_mapping).to_numpy()
        material_feature = np.tan(np.deg2rad(material_feature))

    # Make npz
    trajectories = {}
    if not use_material_feature:
        trajectories[result_dir] = (
            positions_over_time,  # position sequence (timesteps, particles, dims)
            np.full(positions_over_time.shape[1], 6, dtype=int))
    else:
        trajectories[result_dir] = (
            positions_over_time,  # position sequence (timesteps, particles, dims)
            np.full(positions_over_time.shape[1], 6, dtype=int),  # particle type (particles, )
            material_feature)  # particle type (particles, )

    # Create structured array to hold the data
    structured_data = np.empty(len(trajectories), dtype=object)
    for i, value in enumerate(trajectories.values()):
        structured_data[i] = value

    np.savez_compressed(save_path, gns_data=structured_data)
    print(f"Output written to: {save_path}")


convert_hd5_to_npz(
    n_dims=2,
    result_dir='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_0/results/sand3d/',
    mpm_input_path='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_0/mpm.json',
    save_path='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_0/results/debug.npz',
    scale={"origin": [0, 0], "factor": 1/100.0, "timestep_stride": 10},
    use_material_feature=True,
    mpi=32
)



