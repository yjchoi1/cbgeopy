import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import open_npz


def from_h5_to_animation():

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, help="Path to npz file with extension")
    parser.add_argument("--output", type=str, help="Directory to save the animation")
    parser.add_argument("--ndim", type=int, help="dimension of the simulation")
    parser.add_argument("--xboundary", nargs="+", help="x boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument("--yboundary", nargs="+", help="y boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument("--zboundary", nargs="+", default=None, help="z boundary of simulation domain. (e.g., `0.0 1.0`)")
    parser.add_argument('--step_stride', default=5, type=int, help='Stride step for animation frame')
    parser.add_argument('--max_time', default=None, type=int, help='Maximum timestep of h5 file to include in npz file.')
    args = parser.parse_args()

    path = args.path
    output = args.output
    step_stride = args.step_stride
    xboundary = args.xboundary
    yboundary = args.yboundary
    zboundary = args.zboundary
    max_time = args.max_time
    n_dims = args.ndim

    # # Just for debugging
    # path = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/trj.npz"
    # output = "/work2/08264/baagee/frontera/gns-mpm-data/mpm/sand2d_frictions_test/sand2d_inverse_eval31/results/"
    # xboundary = [0, 2]
    # yboundary = [0, 0.75]
    # zboundary = None
    # max_time = None
    # n_dims = 2

    trajectories = open_npz.load_npz_data(path)

    for trj in trajectories:
        positions_over_time = trj[0]

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