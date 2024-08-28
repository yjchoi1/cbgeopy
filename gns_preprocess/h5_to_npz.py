import pathlib
import glob
import re
import json
import h5py
import numpy as np
import os
import pandas as pd
from collections import defaultdict


def concat_h5s_at_t(h5s_at_t, result_dir):
    dfs = []
    for h5_name in h5s_at_t:
        h5_path = os.path.join(result_dir, h5_name)
        df = pd.read_hdf(h5_path, 'table')
        df = df[['id', 'coord_x', 'coord_y', 'coord_z', 'material_id']]
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    return df


def convert_hd5_to_npz(
        n_dims: int,
        result_dir: str,
        mpm_input_path: str,
        save_path: str,
        material_feature: bool,
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
        material_feature ():
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
        time = files_name.split('-')[-1].split('.')[0]
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

    for t, h5s_at_t in files_by_time.items():
        print(t)

        df = concat_h5s_at_t(h5s_at_t, result_dir)

        # number of particles
        n_particles = len(df['coord_x'])

        # positions
        if n_dims == 3:
            positions = df[['coord_x', 'coord_y', 'coord_z']].to_numpy()
        elif n_dims == 2:
            positions = df[['coord_x', 'coord_y']].to_numpy()
        else:
            raise ValueError

        # Check if the particle disappears
        if n_particles_reference is None:
            n_particles_reference = n_particles
        else:
            if n_particles != n_particles_reference:
                raise ValueError(f"Number of particles at timestep {t} is not the same as previous one.")

    # concat positions over time
    positions_over_time = np.concatenate(positions_over_time)

    # material feature
    if material_feature:
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

    # Make npz
    trajectories = {}
    if not material_feature:
        trajectories[result_dir] = (
            positions_over_time,  # position sequence (timesteps, particles, dims)
            np.full(positions_over_time.shape[1], 6, dtype=int))  # particle type (particles, )
    else:
        trajectories[result_dir] = (
            positions_over_time,  # position sequence (timesteps, particles, dims)
            np.full(positions_over_time.shape[1], 6, dtype=int))  # particle type (particles, )


convert_hd5_to_npz(
    n_dims=2,
    result_dir='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_1/results/sand3d/',
    mpm_input_path='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_1/mpm.json',
    save_path='/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_scale_1/results/',
    material_feature=True,
)



