import trimesh
import pandas as pd
import numpy as np
import trimesh
from scipy.spatial import Delaunay, cKDTree
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def trapezoidslope(cell_size, n_particle_per_dim, save_path):

    # Default variables
    particle_distance = cell_size / n_particle_per_dim
    particle_offset = particle_distance / 2
    y_min = 0.0
    y_max = 0.25
    y_interval = particle_distance

    def surface(x):
        if 0 <= x < 0.15:
            z = 0.25
        elif 0.15 <= x < 0.50:
            z = -(5 / 7) * x + 5/14
        else:
            z = 0
        return z

    def surface2(x):
        if 0 <= x < 0.50:
            z = 0.25
        else:
            z = 0
        return z

    def surface3(x):
        if 0 <= x < 0.3:
            z = 0.25
        elif 0.3 <= x < 0.4:
            z = -(5 / 2) * (x - 0.4)
        else:
            z = 0
        return z

    def base(x):
        if 0 < x < 0.05:
            z = 0.25
        elif 0.05 <= x < 0.3:
            z = -x + 0.30
        else:
            z = 0
        return z

    x_base = np.arange(0 + particle_offset, 0.5, particle_distance)
    y_values = np.arange(y_min + particle_offset, y_max, y_interval)

    particles = []  # Coordinate of particles with shape=(n_particles, 3)
    for x in x_base:
        z_min = 0
        z_max = surface3(x)
        zs = np.arange(z_min + particle_offset, z_max, particle_distance)
        for z in zs:
            for y in y_values:
                particles.append((x, y, z))
    particles = np.array(particles)

    pset0 = []
    pset1 = []
    for particle in particles:
        x, y, z = particle[0], particle[1], particle[2]
        if base(x) <= z < surface3(x):
            pset0.append(particle)
        else:
            pset1.append(particle)

    psets = {'0': np.array(pset0), '1': np.array(pset1)}

    # Save the particles to the specified file
    for key, value in psets.items():
        df_particles = pd.DataFrame(value, columns=['x', 'y', 'z'])
        df_particles.to_csv(f'{save_path}/particles_{key}.csv', index=False, header=False)

    return psets


def reduce_mesh(input_path, save_path):
    mesh = trimesh.load_mesh(input_path)

    # Reduce the resolution of the mesh
    simplified_mesh = mesh.simplify_quadric_decimation(face_count=50000)

    # Save the simplified mesh to a new .obj file
    simplified_mesh.export(save_path)


def generate_slope(x, y, slope_angle, z_min, amplitude=1, frequency=1):
    """
    Define a 3D plane with a specified slope angle and sinusoidal variation.

    Parameters:
    - x: array-like, x coordinates
    - y: array-like, y coordinates
    - slope_angle: float, slope angle with respect to the x-axis in degrees
    - amplitude: float, amplitude of the sinusoidal variation
    - frequency: float, frequency of the sinusoidal variation

    Returns:
    - z: array-like, z coordinates
    """
    # Convert slope angle from degrees to radians
    theta = np.radians(slope_angle)

    # Calculate the slope (tangent of the angle)
    slope_x = np.tan(theta)

    # Define the base plane equation z = slope_x * x
    z_base = slope_x * x

    # Add sinusoidal variation
    z_variation = amplitude * np.sin(frequency * x) * np.cos(frequency * y)

    # Combine base plane and variation
    z = z_base + z_variation + z_min

    return z


def gen_surface_mesh(x_range, y_range, grid_distance, height_function, plot=False):
    """

    Args:
        x_range (tuple): range of x values (min, max)
        y_range (tuple): range of y values (min, max)
        grid_distance (float): distance between mesh grid
        height_function (callable): a function of x and y that defines the height at each point
        plot ():

    Returns:
        trimesh.Trimesh object representing the 3D mesh surface
    """

    x = np.arange(x_range[0], x_range[1] + grid_distance, grid_distance)
    y = np.arange(y_range[0], y_range[1] + grid_distance, grid_distance)
    x, y = np.meshgrid(x, y)
    z = height_function(x, y)

    # Flatten the arrays to create points
    points = np.vstack((x.flatten(), y.flatten())).T

    # Delaunay triangulation to create faces
    tri = Delaunay(points)
    faces = tri.simplices

    # Create vertices array
    vertices = np.vstack((x.flatten(), y.flatten(), z.flatten())).T

    # Create the mesh
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    # Plot the mesh using matplotlib
    if plot == True:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles=faces, cmap='viridis', edgecolor='none')
        ax.set_title('3D topography of Underground Layer')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show()

    return mesh


