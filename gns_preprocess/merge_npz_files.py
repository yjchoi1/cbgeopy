import numpy as np
import json
from tqdm import tqdm

# Inputs
bounds = [[0, 300], [0, 150]]
sequence_length = 400
default_connectivity_radius = 2.5
dim = int(2)
material_feature_len = int(1)
dt_mpm = 0.375  # 0.0025
mpm_cell_size = 4  # [0.0125, 0.0125]
nparticles_per_cell = 4  # int(16)
dt_gns = 1.0  # 1.0 is default

result_root_dir = "/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d_layers-random/"  # "./mpm"
result_subdir = "/"
npz_name = "trajectory"  # "mpm-9k-train"
data_tags = [i for i in range(744, 745)] \
            # + [i for i in range(180, 210)] \
            # + [i for i in range(360, 390)] \
            # + [i for i in range(540, 570)] \
            # + [i for i in range(720, 750)]
# excluded_data_tags = [179]
# data_tags = [i for i in data_tags if i not in excluded_data_tags]
save_dir = './'
save_name = "sand2d-layer"

# data containers
trajectories = {}
data_ids = []
# for computing statistics
cumulative_count = 0
cumulative_sum_vel = np.zeros((1, dim))
cumulative_sum_acc = np.zeros((1, dim))
cumulative_sumsq_vel = np.zeros((1, dim))
cumulative_sumsq_acc = np.zeros((1, dim))

# iterate over each simulation
for id in tqdm(data_tags, total=len(data_tags)):

    data_ids.append(id)
    npz_path = f"{result_root_dir}/sim-{id}/{result_subdir}/{npz_name}.npz"  # f"{mpm_dir}/{data_name}/{data_name}.npz"
    data = np.load(npz_path, allow_pickle=True)
    # get trajectory info
    if 'gns_data' in data:
        # for latest npz data
        try:
            data_dict = data['gns_data'].item()
            simulation_id = list(data_dict.keys())[0]
            trajectory = list(data_dict.values())[0]
            trajectories[simulation_id] = trajectory
        except:
            trajectory = data['gns_data'][0]
            trajectories[f"simulation_trajectory_{id}"] = trajectory
    else:
        # for previous npz data
        for simulation_id, trajectory in data.items():  # note, only one trajectory exists, so no need to iterate
            trajectories[f"simulation_trajectory_{id}"] = (trajectory)

    # get positions for each mpm simulation
    positions = trajectory[0]
    if any(dim == 0 for dim in positions.shape):
        raise ValueError(f"The values for the result file are empty in {npz_path}")
    array_shape = positions.shape
    flattened_positions = np.reshape(positions, (-1, array_shape[-1]))

    # compute velocities using finite difference
    # assume velocities before zero are equal to zero
    velocities = np.empty_like(positions)
    velocities[1:] = (positions[1:] - positions[:-1]) / dt_gns
    velocities[0] = 0
    flattened_velocities = np.reshape(velocities, (-1, array_shape[-1]))

    # compute accelerations finite difference
    # assume accelerations before zero are equal to zero
    accelerations = np.empty_like(velocities)
    accelerations[1:] = (velocities[1:] - velocities[:-1]) / dt_gns
    accelerations[0] = 0
    flattened_accelerations = np.reshape(accelerations, (-1, array_shape[-1]))

    # Compute statistics
    cumulative_count += len(flattened_velocities)
    # running sum
    cumulative_sum_vel += np.sum(flattened_velocities, axis=0)
    cumulative_sum_acc += np.sum(flattened_accelerations, axis=0)
    # running sum squared
    cumulative_sumsq_vel += np.sum(flattened_velocities**2, axis=0)
    cumulative_sumsq_acc += np.sum(flattened_accelerations**2, axis=0)
    # statistics for cumulative data
    cumulative_mean_vel = cumulative_sum_vel / cumulative_count
    cumulative_mean_acc = cumulative_sum_acc / cumulative_count
    # cumulative_std_vel = np.sqrt(
    #     (cumulative_sumsq_vel - cumulative_sum_vel ** 2 / cumulative_count) / (cumulative_count - 1))
    # cumulative_std_acc = np.sqrt(
    #     (cumulative_sumsq_acc - cumulative_sum_acc ** 2 / cumulative_count) / (cumulative_count - 1))
    cumulative_std_vel = np.sqrt(
        (cumulative_sumsq_vel/cumulative_count - (cumulative_sum_vel/cumulative_count)**2))
    cumulative_std_acc = np.sqrt(
        (cumulative_sumsq_acc/cumulative_count - (cumulative_sum_acc/cumulative_count)**2))

# Store final statistics
if dim == 2:
    statistics = {
        "mean_velocity_x": float(cumulative_mean_vel[:, 0]),
        "mean_velocity_y": float(cumulative_mean_vel[:, 1]),
        "std_velocity_x": float(cumulative_std_vel[:, 0]),
        "std_velocity_y": float(cumulative_std_vel[:, 1]),
        "mean_accel_x": float(cumulative_mean_acc[:, 0]),
        "mean_accel_y": float(cumulative_mean_acc[:, 1]),
        "std_accel_x": float(cumulative_std_acc[:, 0]),
        "std_accel_y": float(cumulative_std_acc[:, 1])
    }
if dim == 3:
    statistics = {
        "mean_velocity_x": float(cumulative_mean_vel[:, 0]),
        "mean_velocity_y": float(cumulative_mean_vel[:, 1]),
        "mean_velocity_z": float(cumulative_mean_vel[:, 2]),
        "std_velocity_x": float(cumulative_std_vel[:, 0]),
        "std_velocity_y": float(cumulative_std_vel[:, 1]),
        "std_velocity_z": float(cumulative_std_vel[:, 2]),
        "mean_accel_x": float(cumulative_mean_acc[:, 0]),
        "mean_accel_y": float(cumulative_mean_acc[:, 1]),
        "mean_accel_z": float(cumulative_mean_acc[:, 2]),
        "std_accel_x": float(cumulative_std_acc[:, 0]),
        "std_accel_y": float(cumulative_std_acc[:, 1]),
        "std_accel_z": float(cumulative_std_acc[:, 2])
    }

# Print statistics
for key, value in statistics.items():
    print(f"{key}: {value:.7E}")

# Create structured array to hold the data
structured_data = np.empty(len(trajectories), dtype=object)
for i, value in enumerate(trajectories.values()):
    structured_data[i] = value

# Save npz
print(f"Compressing npz...")
np.savez_compressed(f"{save_name}.npz", gns_data=structured_data)
print(f"npz saved at: ./{save_name}.npz")

# Save metadata.json
if dim == 2:
    metadata = {
        "bounds": bounds,
        "sequence_length": sequence_length,
        "default_connectivity_radius": default_connectivity_radius,
        "boundary_augment": 1.0,
        "material_feature_len": material_feature_len,
        "dim": dim,
        "dt": dt_mpm,
        "vel_mean": [statistics["mean_velocity_x"], statistics["mean_velocity_y"]],
        "vel_std": [statistics["std_velocity_x"], statistics["std_velocity_y"]],
        "acc_mean": [statistics["mean_accel_x"], statistics["mean_accel_y"]],
        "acc_std": [statistics["std_accel_x"], statistics["std_accel_y"]],
        "mpm_cell_size": mpm_cell_size,
        "nparticles_per_cell": nparticles_per_cell,
        "data_id": data_ids
    }
elif dim == 3:
    metadata = {
        "bounds": bounds,
        "sequence_length": sequence_length,
        "default_connectivity_radius": default_connectivity_radius,
        "boundary_augment": 1.0,
        "material_feature_len": material_feature_len,
        "dim": dim,
        "dt": dt_mpm,
        "vel_mean": [statistics["mean_velocity_x"], statistics["mean_velocity_y"], statistics["mean_velocity_z"]],
        "vel_std": [statistics["std_velocity_x"], statistics["std_velocity_y"], statistics["std_velocity_z"]],
        "acc_mean": [statistics["mean_accel_x"], statistics["mean_accel_y"], statistics["mean_accel_z"]],
        "acc_std": [statistics["std_accel_x"], statistics["std_accel_y"], statistics["std_accel_z"]],
        "mpm_cell_size": mpm_cell_size,
        "nparticles_per_cell": nparticles_per_cell,
        "data_id": data_ids
    }
else:
    raise ValueError

with open(f"metadata-{save_name}.json", "w") as fp:
    json.dump(metadata, fp, indent=4)
print(f"metadata saved at: ./{save_name}.json")
