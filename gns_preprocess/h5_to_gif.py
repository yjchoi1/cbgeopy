import argparse
import pathlib
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
from collections import defaultdict
import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import utils


def from_h5_to_animation():

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, help="Directory of folder containing mpm h5 files")
    parser.add_argument("--output", type=str, help="Directory to save the animation")
    parser.add_argument("--ndim", type=int, help="dimension of the simulation")
    parser.add_argument("--xboundary", nargs="+", help="x boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument("--yboundary", nargs="+", help="y boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument("--zboundary", nargs="+", default=None, help="z boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument('--step_stride', default=5, type=int, help='Stride step for animation frame')
    parser.add_argument('--max_time', default=None, type=int, help='Maximum timestep of h5 file to include in npz file.')
    args = parser.parse_args()

    result_dir = args.path
    output = args.output
    step_stride = args.step_stride
    xboundary = args.xboundary
    yboundary = args.yboundary
    zboundary = args.zboundary
    max_time = args.max_time
    n_dims = args.ndim

    # # Just for debugging
    # result_dir = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/metadata-sand2d_inverse_eval"
    # output = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/"
    # xboundary = [0, 2]
    # yboundary = [0, 0.75]
    # zboundary = None
    # max_time = None
    # n_dims = 2

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
        positions_list, axis=0)[::step_stride]

    if n_dims == 2:
        # make animation
        fig, ax = plt.subplots()

        def animate(i):
            fig.clear()
            ax = fig.add_subplot(111, aspect='equal', autoscale_on=False)
            ax.set_xlim([float(xboundary[0]), float(xboundary[1])])
            ax.set_ylim([float(yboundary[0]), float(yboundary[1])])
            ax.scatter(positions_over_time[i][:, 0], positions_over_time[i][:, 1], s=1)
            ax.grid(True, which='both')

        ani = animation.FuncAnimation(
            fig, animate, frames=np.arange(0, len(positions_over_time), step_stride), interval=10)

        ani.save(f'{output}/trajectory.gif', dpi=100, fps=30, writer='imagemagick')
        print(f"Animation saved to: {output}")

    if n_dims == 3:
        # make animation
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        def animate(i):
            fig.clear()
            ax = fig.add_subplot(projection='3d', autoscale_on=False)
            ax.set_xlim([float(xboundary[0]), float(xboundary[1])])
            ax.set_ylim([float(yboundary[0]), float(yboundary[1])])
            ax.set_zlim([float(zboundary[0]), float(zboundary[1])])
            ax.scatter(positions_over_time[i][:, 0], positions_over_time[i][:, 1], positions_over_time[i][:, 2], s=1)
            ax.set_box_aspect(
                aspect=(float(xboundary[1]) - float(xboundary[0]),
                        float(yboundary[1]) - float(yboundary[0]),
                        float(zboundary[1]) - float(zboundary[0])))
            ax.view_init(elev=20., azim=i*0.5)
            ax.grid(True, which='both')

        ani = animation.FuncAnimation(
            fig, animate, frames=np.arange(0, len(positions_over_time), step_stride), interval=10)

        ani.save(f'{output}/trajectory.gif', dpi=100, fps=30, writer='imagemagick')
        print(f"Animation saved to: {output}")

if __name__ == "__main__":
    from_h5_to_animation()