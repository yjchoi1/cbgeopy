from typing import List, Dict, Tuple
import pandas as pd
import os


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