import numpy as np
import ezdxf
import trimesh
import alphashape
from shapely.geometry import Point
from scipy.spatial import Delaunay, cKDTree
from shapely.geometry import Polygon
from scipy.interpolate import griddata
import os
import pandas as pd
import demo_utils


def save_script(current_script_path, save_path):
    # Read the content of the specified script
    with open(current_script_path, 'r') as script_file:
        script_content = script_file.read()

    # Save the content to the specified path
    with open(save_path, 'w') as save_file:
        save_file.write(script_content)


def get_z_coordinates(
        mesh: trimesh.Trimesh,
        xy_coords: np.array,
        find_method: str
):
    """
    Get the z-coordinate values on the mesh that correspond to given (x, y) coordinates.

    Parameters:
    - mesh_path: trimesh object
    - xy_coords: np.array of shape (n, 2), array of (x, y) coordinates to get z-coords
    - find_method: 'kdtree', 'linear'

    Returns:
    - np.array of shape (n,), array of z-coordinate values.
    """
    # Extract the vertices of the mesh
    vertices = mesh.vertices

    # Create a KDTree for fast nearest-neighbor lookup
    xy_vertices = vertices[:, :2]
    z_vertices = vertices[:, 2]

    if find_method == "kdtree":

        kdtree = cKDTree(xy_vertices)

        # Find the nearest vertex in the mesh for each (x, y) coordinate
        distances, vertex_indices = kdtree.query(xy_coords)

        # Get the z-values for the corresponding vertices
        z_values = vertices[vertex_indices, 2]

    elif find_method == "linear":
        # Use griddata to interpolate the z-values for the given (x, y) coordinates
        z_values = griddata(
            xy_vertices, z_vertices, xy_coords, method='linear')

    else:
        raise ValueError

    return z_values


def fill_particles_between_mesh(
        lower_mesh: trimesh.Trimesh,
        upper_mesh: trimesh.Trimesh,
        cell_size: list,
        n_particles_per_dim: int,
        z_find_method: str,
        base_find_method: str,
        z_fill_method: str = 'simple',
        initial_alpha: float = 0.1,
        alpha_decay: float = 1
):
    """

    Args:
        base_find_method (str): method to find the base of the area where two surfaces overlap ('alphashape' or 'simple')
        z_find_method (str): method to find z-coordinate of mesh
        z_fill_method (str): method to fill between lower and upper z-coordinates
        lower_mesh (trimesh.Trimesh):
        upper_mesh (trimesh.Trimesh):
        cell_size (list):
        initial_alpha (float): = 0.1
        alpha_decay (float):
        n_particles_per_dim (int):

    Returns:

    """
    # Error handling
    # if not all(i == cell_size[0] for i in cell_size):
    #     raise NotImplementedError("Current version only support the same element size for all dimension")

    # Define spacings
    particle_distance = cell_size[0] / n_particles_per_dim
    particle_offset_distance = particle_distance / 2

    # Project the topography mesh to xy plain
    projected_vertices = upper_mesh.vertices[:, :2]

    # Define candidate base points
    # define min anx max grid coord
    # min_grid_values = np.floor(np.min(projected_vertices, axis=0) / cell_size[0]) * cell_size[0]
    # max_grid_values = np.ceil(np.max(projected_vertices, axis=0) / cell_size[0]) * cell_size[0]

    # define min and max particle coord
    min_particle_values = np.min(projected_vertices, axis=0) + particle_offset_distance
    max_particle_values = np.max(projected_vertices, axis=0) - particle_offset_distance

    # define range
    candidate_point_range = list(zip(min_particle_values, max_particle_values))
    rounded_array_processed = []
    for x, y in candidate_point_range:
        x_rounded = round((x - particle_offset_distance) / particle_distance) * particle_distance + particle_offset_distance
        y_rounded = round((y - particle_offset_distance) / particle_distance) * particle_distance + particle_offset_distance
        rounded_array_processed.append((x_rounded, y_rounded))

    # Generate candidate base particle points
    candidate_points = generate_points(
        rounded_array_processed, distance=particle_distance)

    if base_find_method == "alphashape":
        # TODO: add visualization for perimeter check
        # Create a 2D polygon encompassing the projected vertices
        max_iterations = 100
        for i in range(max_iterations):
            alpha = initial_alpha / (alpha_decay * (i + 1))
            alpha_shape = alphashape.alphashape(projected_vertices, alpha=alpha)
            if isinstance(alpha_shape, Polygon):
                break
        else:
            raise ValueError("Could not find a Polygon for projected base in alpha shape")

        # Check which grid points fall inside the projected polygon
        inside_check = [alpha_shape.contains(Point(point)) for point in candidate_points]
        base_particles = candidate_points[inside_check]

    elif base_find_method == "simple":
        base_particles = candidate_points

    else:
        raise ValueError(f"Check `base_find_method`")

    # Define lower and upper z bounds and populate particles inbetween.
    lower_z = get_z_coordinates(
        lower_mesh, base_particles, z_find_method)
    upper_z = get_z_coordinates(
        upper_mesh, base_particles, z_find_method)
    criteria = upper_z - lower_z
    lower_z_processed = np.where(criteria <= 0, 0, lower_z)
    upper_z_processed = np.where(criteria <= 0, 0, upper_z)
    z_ranges = np.column_stack(
        (lower_z_processed, upper_z_processed))
    material_points = fill_particles_inbetween(
        base_particles, z_ranges, particle_distance, method=z_fill_method)

    # ----------------------------------------
    # # Populate particles between two surfaces
    # lower_z = get_z_coordinates(
    #     lower_mesh, base_particles, z_find_method)
    # lower_zgrid = np.floor(lower_z / cell_size[0]) * cell_size[0]
    # lower_zparticle = lower_zgrid + particle_offset_distance
    #
    # upper_z = get_z_coordinates(
    #     upper_mesh, base_particles, z_find_method)
    # upper_zgrid = np.round(upper_z / cell_size[0]) * cell_size[0]
    # upper_zparticle = upper_zgrid - particle_offset_distance
    #
    # criteria = upper_zgrid - lower_zgrid
    #
    # lower_zgrid_processed = np.where(criteria <= 0, 0, lower_zgrid)
    # upper_zgrid_processed = np.where(criteria <= 0, 0, upper_zgrid)
    # zgrid_ranges = np.column_stack(
    #     (lower_zgrid_processed, upper_zgrid_processed))
    #
    # material_points = fill_particles_inbetween(
    #     base_particles, zgrid_ranges, particle_distance)
    # ----------------------------------------

    return material_points


def fill_particles_inbetween(
        xy_coords: np.array,
        z_ranges: np.array,
        particle_distance: float,
        method: str = "simple"
):
    """

    Args:
        xy_coords (np.array): shape=(n_coords, 2)
        z_ranges (np.array): shape=(n_xy_coords, 2) where z_ranges[:, 0] is the lower bound & z_ranges[:, 1] is upper bound
        cell_size (list): [x_len, y_len, z_len]
        particle_distance (float): default distance between particles
        method (str): method to fill particles between z_ranges: "simple", "round"
            If `simple`, it uses the exact z-coordinate of mesh.
            If `round`, it uses the nearest particle grid points for particle generation.

    Returns:

    """
    # particle_distance_offset
    particle_offset_distance = particle_distance / 2

    # Initialize an empty list to store particles
    particles_list = []

    for i in range(len(xy_coords)):
        if method == "simple":
            lower_z = z_ranges[i, 0]
            upper_z = z_ranges[i, 1]
            inside_z_points = np.arange(lower_z + particle_offset_distance, upper_z, particle_distance)
        elif method == "round":
            lower_z = np.round(
                (z_ranges[i, 0] - particle_offset_distance) / particle_distance) * particle_distance + particle_offset_distance
            upper_z = np.round(
                (z_ranges[i, 1] - particle_offset_distance) / particle_distance) * particle_distance + particle_offset_distance
            inside_z_points = np.arange(lower_z, upper_z, particle_distance)
        else:
            raise ValueError(f"Method `{method}` is not valid option")


        num_z_points = len(inside_z_points)

        # Repeat the xy coordinates for the number of z values
        xy_repeated = np.repeat(xy_coords[i].reshape(1, 2), num_z_points, axis=0)

        # Stack the xy coordinates and z values
        particles = np.column_stack((xy_repeated, inside_z_points))

        # Append to the list
        particles_list.append(particles)

        # TODO: potential improvements
        # - If `lower_z` is 0, don't concatenate
        # - Only concatenate `upper_z` if the distance between the last element of `main_particles` exceeds a certain threshold

        # ----------------------------------------
        # # Create z values for the current point
        # z_values = np.arange(
        #     z_ranges[i, 0] + distance/2,
        #     z_ranges[i, 1] - distance/2,
        #     distance)
        # num_z_values = len(z_values)
        #
        # # Repeat the xy coordinates for the number of z values
        # xy_repeated = np.repeat(xy_coords[i].reshape(1, 2), num_z_values, axis=0)
        #
        # # Stack the xy coordinates and z values
        # particles = np.column_stack((xy_repeated, z_values))
        #
        # # Append to the list
        # particles_list.append(particles)
        # ----------------------------------------

    # Concatenate all particles into a single array
    all_particles = np.vstack(particles_list)

    return all_particles


def generate_points(
        ranges: list,
        distance: float):
    """

    Args:
        ranges (list): [[x_min, x_max], [y_min, y_max] [z_min, z_max]]
        distance (float):
        round_unit (float):

    Returns:
        coordinate of points with shape=(n_points, n_dims)

    """
    n_dims = len(ranges)
    if n_dims not in [2, 3]:
        raise ValueError("The function only supports 2D and 3D ranges.")

    coordinates = []
    for r in ranges:
        min_val, max_val = r
        coords = np.arange(min_val, max_val + distance, distance)
        coordinates.append(coords)

    points = np.array(np.meshgrid(*coordinates)).T.reshape(-1, n_dims)

    return points


def obj2mesh(file_path):
    vertices = []
    faces = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'v':
                # Convert vertex coordinates to float and add to the list
                vertex = list(map(float, parts[1:]))
                vertices.append(vertex)
            elif parts[0] == 'f':
                # Convert face indices to integer and correct for 0-based indexing
                face = [int(index) - 1 for index in parts[1:]]
                faces.append(face)

    # Convert lists to numpy arrays
    vertices_array = np.array(vertices)
    faces_array = np.array(faces, dtype=int)

    return vertices_array, faces_array


# Extracting points from various entities in the DXF file
def dxf2points(dxf_file):
    # Read the DXF file
    dxf_file = 'topo_surface_Fundao.dxf'
    doc = ezdxf.readfile(dxf_file)

    points = []
    # Checking for different entity types and extracting points
    for entity in doc.modelspace().query('POLYLINE VERTEX POINT 3DFACE'):
        if entity.dxftype() == 'POLYLINE':
            # print("polyline")
            for vertex in entity.vertices:
                points.append([vertex.dxf.location.x, vertex.dxf.location.y, vertex.dxf.location.z])
        elif entity.dxftype() == 'VERTEX':
            # print("vertex")
            points.append([entity.dxf.location.x, entity.dxf.location.y, entity.dxf.location.z])
        elif entity.dxftype() == 'POINT':
            # print("point")
            points.append([entity.dxf.location.x, entity.dxf.location.y, entity.dxf.location.z])
        elif entity.dxftype() == '3DFACE':
            # print("3dface")
            points.append([entity.dxf.vtx0.x, entity.dxf.vtx0.y, entity.dxf.vtx0.z])
            points.append([entity.dxf.vtx1.x, entity.dxf.vtx1.y, entity.dxf.vtx1.z])
            points.append([entity.dxf.vtx2.x, entity.dxf.vtx2.y, entity.dxf.vtx2.z])
            points.append([entity.dxf.vtx3.x, entity.dxf.vtx3.y, entity.dxf.vtx3.z])
    return np.array(points)


def get_h5(uuid_dir, timestep, n_mpis):
    # Create an empty list to store DataFrames
    dfs = []

    # Iterate over different files from MPI and append to list
    for i in range(n_mpis):
        file = f'particles-{i}_{n_mpis}-{timestep}.h5'  # ex) particles-26_32-0120000.h5
        h5_path = os.path.join(uuid_dir, file)
        dfs.append(pd.read_hdf(h5_path, 'table'))

        # Concatenate all DataFrames
        df = pd.concat(dfs, ignore_index=True)

    return df
