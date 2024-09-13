import glob
import json
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple
import argparse
import utils


LINEAR_ELASTIC_FEATURE = 1.0
KINEMATIC_PARTICLE = 6
NON_KINEMATIC_PARTICLE = 3


def convert_hd5_to_npz(
        n_dims: int,
        result_dir: str,
        mpm_input_path: str,
        save_path: str,
        use_material_feature: bool,
        scale: Dict,
        max_time: int = None
):
    """
    Read the results of CB-Geo MPM and convert it to npz format, which is compatible with GNS.
    It can
    Args:
        n_dims (int): dimensionality of the simulation
        result_dir (str): path to directory where h5 files are located
        mpm_input_path (str): path to mpm.json file with extension
        save_path (str): path to save the result npz file with extension
        use_material_feature (bool): whether to add material feature in npz file
        scale (Dict): You can scale the positions of particles based on a specified origin with scale factor
            Format is dict: {"origin": [x, y, z], "factor": 2.0}
        max_time (int): Maximum timestep of h5 file to include in npz file.

    Returns:
        `.npz`

    """
    print(f"Process {result_dir}...")
    file_paths = glob.glob(f'{result_dir}/particles*.h5')
    files_names = [file_path.split("/")[-1] for file_path in file_paths]

    # Initialize a dictionary to hold the groups
    files_by_time = defaultdict(list)

    # Group the file names by the same timestep
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
            utils.process_time_step, t, h5s_at_t, result_dir, n_dims)
            for t, h5s_at_t in files_by_time.items()]
        for future in as_completed(futures):
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
    print(f"Sampled {len(positions_over_time)} position sequence...")

    positions_over_time = utils.scale_positions(
        positions_over_time, scale["origin"], scale["factor"])

    # fig, ax = plt.subplots()
    # plottime = -1
    # ax.scatter(positions_over_time[plottime, :, 0], positions_over_time[plottime, :, 1])
    # plt.show()

    # material feature
    if use_material_feature:
        # open mpm.json
        with open(mpm_input_path, 'r') as file:
            mpm_json = json.load(file)

        # get material id and associated properties
        materials = mpm_json["materials"]

        # get the map between material_id and corresponding property
        particle_material_mapping = {}  # {material_id: phi_deg}
        particle_type_mapping = {}  # {material_id: material_type}
        for material in materials:
            if "MohrCoulomb" in material["type"] and "friction" in material:
                particle_material_mapping[material["id"]] = material["friction"]
                particle_type_mapping[material["id"]] = KINEMATIC_PARTICLE
            elif "LinearElastic" in material["type"]:
                # Set an arbitrary feature for the LE model. Here, we use 1.0.
                particle_material_mapping[material["id"]] = LINEAR_ELASTIC_FEATURE
                particle_type_mapping[material["id"]] = NON_KINEMATIC_PARTICLE
            else:
                raise NotImplemented("Not supported material")

        # make material_feature based on the id.
        friction_angle_deg = df['material_id'].map(particle_material_mapping).to_numpy()
        material_feature = np.tan(np.deg2rad(friction_angle_deg))

        # make particle_type feature
        particle_types = df['material_id'].map(particle_type_mapping).to_numpy()

    # Make npz
    trajectories = {}
    if not use_material_feature:
        trajectories[result_dir] = (
            positions_over_time,  # position sequence (timesteps, particles, dims)
            np.full(positions_over_time.shape[1], 6, dtype=int))
    else:
        trajectories[result_dir] = (
            positions_over_time.astype("float32"),  # position sequence (timesteps, particles, dims)
            particle_types.astype("int32"),  # particle type (particles, )
            material_feature.astype("float32"))  # particle type (particles, )

    # Create structured array to hold the data
    structured_data = np.empty(len(trajectories), dtype=object)
    for i, value in enumerate(trajectories.values()):
        structured_data[i] = value

    np.savez_compressed(save_path, gns_data=structured_data)
    print(f"Output written to: {save_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert mpm result in hdf5 trajectories to npz.')
    parser.add_argument(
        '--n_dims', type=int,
        help="MPM simulation dimensions")
    parser.add_argument(
        '--result_dir', type=str,
        help="Directory where h5 files are located.")
    parser.add_argument(
        '--mpm_input_path', type=str,
        help="Path to mpm.json file. This is used to get material feature")
    parser.add_argument(
        '--save_path', type=str,
        help="Path to save npz file with extension `.npz")
    parser.add_argument(
        '--scale_origin', type=float, nargs=2, required=False,
        help="Origin values where the scaling is based on. Default is [0.0, 0.0, (0.0)]"
             "Pass as two space-separated values (e.g., --origin 0.0 0.0).")
    parser.add_argument(
        '--scale_factor', type=float, required=False,
        help="Factor value for the scaling. Default is 1.0."
             " (e.g., --factor 0.01 for 1/100.0).")
    parser.add_argument(
        '--timestep_stride', type=int, required=False,
        help="Timestep stride value. Default is 1.0."
             " (e.g., --timestep_stride 10).")
    args = parser.parse_args()

    # Just for debug
    # n_dims = 2
    # result_dir = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/metadata-sand2d_inverse_eval"
    # mpm_input_path = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/mpm_input.json"
    # save_path = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/trj.npz"

    scale = {
        "origin": [0, 0] if args.scale_origin is None else args.scale_origin,
        "factor": [1.0, 1.0] if args.scale_factor is None else args.scale_factor,
        "timestep_stride": int(1) if args.timestep_stride is None else args.timestep_stride
    }

    convert_hd5_to_npz(
        n_dims=args.n_dims,
        result_dir=args.result_dir,
        mpm_input_path=args.mpm_input_path,
        save_path=args.save_path,
        scale=scale,
        use_material_feature=True)



