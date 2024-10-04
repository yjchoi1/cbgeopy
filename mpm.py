import numpy as np
import pandas
import plotly.graph_objects as go
import pandas as pd
import os
import json
import trimesh
import utils
import random
import argparse
from scipy.spatial import cKDTree, KDTree
from scipy import interpolate
from typing import List, Dict, Callable


# Define constants for axes
AXES_3D = ["x", "y", "z"]
AXES_2D = ["x", "y"]


class MPMConfig:
    def __init__(
            self,
            domain_origin,
            domain_length,
            title="mpm_input_config"
    ):
        # Mesh
        self.initial_stresses = None
        self.materials = None
        self.cell_size = None
        self.n_cells_per_dim = None
        self.nnode = None
        self.nele = None
        self.mesh_info = None
        self.mesh_coord_base = None

        # Entity
        self.entity_sets = {}

        # Particle
        self.particle_group_id = 0
        self.particles_count = 0
        self.particle_groups = {}

        # General
        self.direction_mapping = {"x": 0, "y": 1, "z": 2}
        self.ndims = len(domain_length)
        self.domain_origin = domain_origin
        self.domain_length = domain_length
        self.domain_ranges = [
            (origin, origin + length) for origin, length in zip(domain_origin, domain_length)]

        # Input json file
        self.mpm_json = {
            "title": title,
            "mesh": {
                "mesh": "mesh.txt",
                "entity_sets": "entity_sets.json",
                "boundary_conditions": {},
                "cell_type": "ED3H8" if self.ndims == 3 else "ED2Q4",
                "isoparametric": False,
                "check_duplicates": True,
                "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                "node_type": "N3D" if self.ndims == 3 else "N2D"
            },
            "particles": []
        }

    def add_mesh(
            self,
            n_cells_per_dim: List[int],
            outer_cell_thickness: float = 0,
    ):
        """
        Make mesh coordinate array & cell group
        Args:
            outer_cell_thickness (float): If provided, the outer mesh will be added to the mesh defined by
                domain_ranges. Therefore, The actual domain is
                [domain_ranges[0] - outer_cell_thickness, domain_ranges[1] + outer_cell_thickness]
                If not provided, the actual domain is
                [domain_ranges[0], domain_ranges[1]]
            n_cells_per_dim (list): [nx, ny, nz]

        Returns:
            None

        """
        self.n_cells_per_dim = n_cells_per_dim

        # # Error handling
        # if not all(i == n_cells_per_dim[0] for i in n_cells_per_dim):
        #     raise NotImplementedError("Current version only support the same element size for all dimension")

        # Make coordinate base for each axis
        self.mesh_coord_base = []
        self.cell_size = []
        for dim, domain_range in enumerate(self.domain_ranges):
            if outer_cell_thickness == 0:
                coord_base = np.linspace(
                    domain_range[0], domain_range[1], n_cells_per_dim[dim] + 1, endpoint=True)
            elif outer_cell_thickness > 0:
                first_cell = np.array([domain_range[0] - outer_cell_thickness])
                end_cell = np.array([domain_range[1] + outer_cell_thickness])
                coord_base = np.concatenate((
                    first_cell,
                    np.linspace(domain_range[0], domain_range[1], n_cells_per_dim[dim]+1, endpoint=True),
                    end_cell))
            else:
                raise ValueError("Outer cell thickness cannot be negative value")
            self.mesh_coord_base.append(coord_base)

            # Calculate the cell size for the current dimension
            cell_size_per_dim = (coord_base.max() - coord_base.min()) / n_cells_per_dim[dim]
            self.cell_size.append(cell_size_per_dim)

        # Create node coordinates
        if self.ndims == 3:
            node_coords = []
            for z in self.mesh_coord_base[2]:
                for y in self.mesh_coord_base[1]:
                    for x in self.mesh_coord_base[0]:
                        node_coords.append([x, y, z])
            node_coords = np.array(node_coords)

            # Compute the number of nodes and elements (cells) for each dimension
            # - node
            nnode_x = len(self.mesh_coord_base[0])
            nnode_y = len(self.mesh_coord_base[1])
            nnode_z = len(self.mesh_coord_base[2])
            self.nnode = nnode_x * nnode_y * nnode_z
            # - element
            nele_x = nnode_x - 1
            nele_y = nnode_y - 1
            nele_z = nnode_z - 1
            self.nele = nele_x * nele_y * nele_z
            nnode_in_ele = 8  # hardcoded for cube for 3-D case

            # Make cell groups
            cells = np.empty((int(self.nele), int(nnode_in_ele)))
            i = 0
            for elz in range(int(nele_z)):
                for ely in range(int(nele_y)):
                    for elx in range(int(nele_x)):
                        # cell index starts from 1 not 0, so there is "1+" at first
                        cells[i, 0] = nnode_x * nnode_y * elz + ely * nnode_x + elx
                        cells[i, 1] = nnode_x * nnode_y * elz + ely * nnode_x + 1 + elx
                        cells[i, 2] = nnode_x * nnode_y * elz + (ely + 1) * nnode_x + 1 + elx
                        cells[i, 3] = nnode_x * nnode_y * elz + (ely + 1) * nnode_x + elx
                        cells[i, 4] = nnode_x * nnode_y * (elz + 1) + ely * nnode_x + elx
                        cells[i, 5] = nnode_x * nnode_y * (elz + 1) + ely * nnode_x + 1 + elx
                        cells[i, 6] = nnode_x * nnode_y * (elz + 1) + (ely + 1) * nnode_x + 1 + elx
                        cells[i, 7] = nnode_x * nnode_y * (elz + 1) + (ely + 1) * nnode_x + elx
                        i += 1
            cells = cells.astype(int)

        else:  # for 2d
            # Create node coordinates
            node_coords = []
            for y in self.mesh_coord_base[1]:
                for x in self.mesh_coord_base[0]:
                    node_coords.append([x, y])
            node_coords = np.array(node_coords)

            # Compute the number of nodes and elements for each dimension
            # - node
            nnode_x = len(self.mesh_coord_base[0])
            nnode_y = len(self.mesh_coord_base[1])
            nnode_z = None
            self.nnode = nnode_x * nnode_y
            # - element
            nele_x = nnode_x - 1
            nele_y = nnode_y - 1
            self.nele = nele_x * nele_y
            nnode_in_ele = 4  # hardcoded for cube for 3-D case

            # Define cell groups consisting of four nodes
            cells = np.empty((int(self.nele), int(nnode_in_ele)))
            i = 0
            for ely in range(int(nele_y)):
                for elx in range(int(nele_x)):
                    # cell index starts from 1 not 0, so there is "1+" at first
                    cells[i, 0] = nnode_x * ely + elx
                    cells[i, 1] = nnode_x * ely + elx + 1
                    cells[i, 2] = nnode_x * (ely + 1) + elx + 1
                    cells[i, 3] = nnode_x * (ely + 1) + elx
                    i += 1
            cells = cells.astype(int)

        # Aggregate data in dict
        self.mesh_info = {
            "node_coords": node_coords,
            "n_node_x": nnode_x,
            "n_node_y": nnode_y,
            "n_node_z": nnode_z,
            "cell_groups": cells
        }

    def add_particles_csv(
            self,
            path,
            material_id,
            particle_group_id=None):
        """

        Args:
            material_id (int): material id associated with this particle group
            particle_group_id (int): particle group id to be associated with this particles
            path (str): csv file path

        Returns:

        """
        # Assign a particle group id
        if particle_group_id is None:
            self.particle_group_id += 1
        else:
            self.particle_group_id = particle_group_id

        # Read the particle CSV file into a pandas DataFrame
        df = pd.read_csv(path)
        # Convert the DataFrame to a NumPy array
        particles = df.to_numpy()

        # Store
        self.particle_groups[self.particle_group_id] = {}
        self.particle_groups[self.particle_group_id]['particles'] = particles
        self.particle_groups[self.particle_group_id]['id'] = list(
            range(self.particles_count, self.particles_count + len(particles)))

        # Update current particle count
        self.particles_count += len(particles)

        # Set config
        self.mpm_json["particles"].append(
            {
                "generator": {
                    "check_duplicates": True,
                    "location": f"particles_{self.particle_group_id}.txt",
                    "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                    "pset_id": self.particle_group_id,
                    "particle_type": "P3D" if self.ndims == 3 else "P2D",
                    "material_id": material_id,
                    "type": "file"
                }
            }
        )

    def add_particles_from_lines(
            self,
            layer_info: List,
            n_particle_per_cell: int,
            randomness: float = None
    ):
        """

        Args:
            layer_info (List): a list of dict {"line_points": [a list of points], "material_id" int}
                * "line_points" contains the points that comprise a line that defines the upper boundary of layer.
                    It should be defined for the entire x-domain range.
                    For example, [[[0, 0], [1.0, 0]], [[0, 0], [0.3, 0.3], [0.7, 0.1], [1.0, 0]], [[0, 0.5], [0.1]]]
                * "material_id" is the material id associated with this layer.
                * "particle_group_id": particle_group_id (int): particle group id to be associated with this particles
            n_particle_per_cell (int): number of particles per cell per dimension
            randomness ():

        Returns:

        """
        if self.ndims == 3:
            raise ValueError("This feature is only for 2D domain")

        # Particle config for whole domain
        particle_distance = self.cell_size[0] / n_particle_per_cell
        particle_offset_distance = particle_distance / 2
        particle_ranges = [
            (origin + particle_offset_distance, origin + length - particle_offset_distance)
            for origin, length in zip(self.domain_origin, self.domain_length)]

        # Create particle range arrays that cover the whole domain
        x_coords = np.arange(
            particle_ranges[0][0], particle_ranges[0][1] + particle_offset_distance, particle_distance)
        y_coords = np.arange(
            particle_ranges[1][0], particle_ranges[1][1] + particle_offset_distance, particle_distance)

        # Generate the candidate particles
        xx, yy = np.meshgrid(x_coords, y_coords)
        candidate_particles = np.vstack((xx.ravel(), yy.ravel())).T

        # Initialize lower boundary (for the first layer it is assumed to be y=0)
        lower_boundary_y = np.zeros_like(candidate_particles[:, 0])

        # Iterate over layers to define upper boundaries
        for i, layer in enumerate(layer_info):
            # Assign a particle group id
            if "particle_group_id" in layer:
                self.particle_group_id = layer["particle_group_id"]
            else:
                self.particle_group_id += 1

            # Create line function for the upper boundary of the layer
            points = np.array(layer["line_points"])
            x_points, y_points = points[:, 0], points[:, 1]
            upper_boundary_func = interpolate.interp1d(
                x_points, y_points, kind='linear', fill_value="extrapolate")

            # Determine the y-values of the upper boundary at all x-coordinates
            upper_boundary_y = upper_boundary_func(candidate_particles[:, 0])

            # Find particles that are between the lower and upper boundaries
            mask = (candidate_particles[:, 1] < upper_boundary_y) & (candidate_particles[:, 1] >= lower_boundary_y)
            particles = candidate_particles[mask]

            # Disturb particles
            if randomness is not None:
                particles += np.random.uniform(
                    -particle_offset_distance * randomness,
                    particle_offset_distance * randomness,
                    particles.shape)

            # the current upper boundary becomes lower boundary of the next iteration
            lower_boundary_y = upper_boundary_y

            # Store
            self.particle_groups[self.particle_group_id] = {}
            self.particle_groups[self.particle_group_id]['particles'] = particles
            self.particle_groups[self.particle_group_id]['id'] = list(
                range(self.particles_count, self.particles_count + len(particles)))

            # Update current particle count
            self.particles_count += len(particles)

            # Set config
            self.mpm_json["particles"].append(
                {
                    "generator": {
                        "check_duplicates": True,
                        "location": f"particles_{self.particle_group_id}.txt",
                        "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                        "pset_id": self.particle_group_id,
                        "particle_type": "P3D" if self.ndims == 3 else "P2D",
                        "material_id": layer["material_id"],
                        "type": "file"
                    }
                }
            )



    def add_particles_from_topography(
            self,
            lower_topography: trimesh.Trimesh,
            upper_topography: trimesh.Trimesh,
            n_particle_per_cell: int,
            material_id: int,
            z_find_method: str,
            base_find_method: str,
            z_fill_method: str = 'simple',
            randomness: float = None,
            overlap_tolerance: float = None,
            particle_group_id: int = None
    ):
        """

        Args:
            randomness (float): disturb particles with uniform randomness factor
            overlap_tolerance (float):
            lower_topography (trimesh.Trimesh): mesh that defines the surface of the upper layer topography
            upper_topography (trimesh.Trimesh):  mesh that defines the surface of the upper layer topography
            n_particle_per_cell (int): number of particles per cell per dimension
            material_id (int): id of material that you want to assign to this particle set
            z_find_method (str): method to find z-coordinate of mesh ('' or '')
            base_find_method (str): method to find the base of the area where two surfaces overlap ('alphashape' or 'simple')
                If 'alphashape', it defines base by trying to estimate the x-y plane projection area using alphashape.
                if 'simple, it defines base using max and min x, y values. The base will be restricted to square shape.
            z_fill_method (str): method to fill between lower and upper z-coordinates ('simple' or 'round')
                If 'simple', it uses the exact z-coordinate of mesh.
                If 'round', it uses the nearest particle grid points for particle generation.
            particle_group_id (int): particle group id to be associated with this particles

        Returns:

        """

        if self.ndims == 2:
            raise ValueError("This feature is only for 3D domain")

        # Particle config
        particle_distance = self.cell_size[0] / n_particle_per_cell
        particle_offset_distance = particle_distance / 2

        # Assign a particle group id
        if particle_group_id is None:
            self.particle_group_id += 1
        else:
            self.particle_group_id = particle_group_id

        # Fill particles between two meshes
        new_particles = utils.fill_particles_between_mesh(
            lower_topography,
            upper_topography,
            self.cell_size,
            n_particle_per_cell,
            z_find_method,
            base_find_method,
            z_fill_method
        )

        # Disturb particles
        if randomness is not None:
            new_particles += np.random.uniform(
                -particle_offset_distance * randomness,
                particle_offset_distance * randomness,
                new_particles.shape)

        # overlap check with previous particles
        overlap_count = 0

        if len(self.particle_groups) > 0 and overlap_tolerance is not None:
            print(f"Overlap checking for particle group {particle_group_id}")
            # Collect all existing particles into a single array for KDTree
            all_existing_particles = np.vstack(
                [group['particles'] for group in self.particle_groups.values()]
            )

            # Build KDTree with existing particles
            tree = KDTree(all_existing_particles)

            # Check for overlapping particles
            non_overlapping_particles = []
            for particle in new_particles:
                dist, _ = tree.query(particle, distance_upper_bound=overlap_tolerance)
                if dist == float('inf'):  # If distance is inf, it means no close neighbor was found
                    non_overlapping_particles.append(particle)
                else:
                    overlap_count += 1

            non_overlapping_particles = np.array(non_overlapping_particles)

        else:
            # If only one particle group exists, all new particles are non-overlapping
            non_overlapping_particles = new_particles

        # Store
        self.particle_groups[self.particle_group_id] = {}
        self.particle_groups[self.particle_group_id]['particles'] = non_overlapping_particles
        self.particle_groups[self.particle_group_id]['id'] = list(
            range(self.particles_count, self.particles_count + len(non_overlapping_particles))
        )

        # Update current particle count
        self.particles_count += len(non_overlapping_particles)
        print(f'Number of overlapping particles found: {overlap_count}')

        # Set config
        self.mpm_json["particles"].append(
            {
                "generator": {
                    "check_duplicates": True,
                    "location": f"particles_{self.particle_group_id}.txt",
                    "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                    "pset_id": self.particle_group_id,
                    "particle_type": "P3D" if self.ndims == 3 else "P2D",
                    "material_id": material_id,
                    "type": "file"
                }
            }
        )

    def define_particle_entity(self):
        """
        Write the indices of particle set & indices of particles for current particle set.
        It is used for mpm input to define different materials for each particle set in mpm solver.
        """
        self.entity_sets['particle_sets'] = []
        for set_id, particle_dict in self.particle_groups.items():
            self.entity_sets["particle_sets"].append({
                "id": set_id,  # index of particle set
                "set": particle_dict['id'],  # index of particles for current set
            })

    def remove_overlapping_particles(self, overlap_tolerance):
        """
        Iterate over all particles in particle groups and remove the overlapping particles
            when the particles in the current group overlaps the previous cumulative particles
        This reorders the particle ids, i.e., `particle_groups['particle_group_id']['id']`
        Args:
            overlap_tolerance ():

        Returns:

        """

        if len(self.particle_groups) == 0:
            raise ValueError("There are no existing particle groups")

        overlap_count = 0
        self.particles_count = 0
        particle_group_count = 0
        current_existing_particles = None

        for pset_id, pdata in self.particle_groups.items():
            new_particles = pdata['particles']

            # For the first particle group, all new particles are non-overlapping
            if particle_group_count == 0:
                non_overlapping_particles = new_particles
                current_existing_particles = new_particles

            # overlap check with previous particles
            else:
                print(f"Overlap checking for particle group {pset_id}")

                # Build KDTree with existing particles
                tree = KDTree(current_existing_particles)

                # Check for overlapping particles
                non_overlapping_particles = []
                for particle in new_particles:
                    dist, _ = tree.query(particle, distance_upper_bound=overlap_tolerance)
                    if dist == float('inf'):  # If distance is inf, it means no close neighbor was found
                        non_overlapping_particles.append(particle)
                    else:
                        overlap_count += 1

                non_overlapping_particles = np.array(non_overlapping_particles)

                # Update existing particles with non-overlapping ones
                if len(non_overlapping_particles) > 0:
                    current_existing_particles = np.vstack((current_existing_particles, non_overlapping_particles))

            particle_group_count += 1

            # Update particle count after filtering
            self.particles_count += len(non_overlapping_particles)

            # Store
            self.particle_groups[pset_id]['particles'] = non_overlapping_particles
            self.particle_groups[pset_id]['id'] = list(
                range(self.particles_count, self.particles_count + len(non_overlapping_particles))
            )

            print(f'Number of overlapping particles found: {overlap_count}')

    def add_cell_entity(
            self,
            nset_id: int,
            ranges: list):
        """

        Args:
            ranges (list): [[x_min, x_max], [y_min, y_max], [z_min, z_max]].

        Returns:

        """
        # TODO: set id should be properly considered. Currently, it is first auto created by boundary entity
        if nset_id <= 5:
            raise ValueError("`nset_id` should be larger than 5 since 0 to 5 is already occupied by the boundary node sets")

        node_coords = self.mesh_info['node_coords']
        nodes = []
        coords = []

        # Ensure ranges are properly formatted according to dimensionality
        if len(ranges) != self.ndims:
            raise ValueError(f"Expected {self.ndims} ranges, but got {len(ranges)}")

        # Determine nodes within the specified ranges
        # ranges[dim][0] and ranges[dim][1] is the min and max range in dim.
        for i, point in enumerate(node_coords):
            if all(ranges[dim][0] <= point[dim] <= ranges[dim][1] for dim in range(self.ndims)):
                nodes.append(i)
                coords.append(point)

        # Add the node set to the entity sets
        self.entity_sets["node_sets"].append({
            "id": nset_id,  # index of node set
            "set": nodes,  # index of nodes
        })

    def add_particles_cube(
            self,
            cube_origin,
            cube_length,
            material_id,
            n_particle_per_cell,
            randomness: float = None,
            particle_group_id=None):
        """
        Adds a group of particles within a defined cubic region of the simulation domain.

        Args:
            randomness (float): disturb particles with uniform randomness factor
            cube_origin (List[float]): The lower edge of the cube, [x_min, y_min, (z_min)].
            cube_length (List[float]): The lengths of the cube along each axis, [x_len, y_len, (z_len)].
            material_id (int): The material ID associated with this particle group.
            n_particle_per_cell (int): Number of particles per dimension in each cell.
            particle_group_id (Optional[int]): Particle group ID to be associated with these particles. Auto-increments if None.

        Raises:
            ValueError: If an unsupported number of dimensions is provided.
        """

        # Assign a particle group id
        if particle_group_id is None:
            self.particle_group_id += 1
        else:
            self.particle_group_id = particle_group_id

        # Particle config
        particle_distance = self.cell_size[0] / n_particle_per_cell
        particle_offset_distance = particle_distance / 2
        particle_ranges = [
            (origin + particle_offset_distance, origin + length - particle_offset_distance)
            for origin, length in zip(cube_origin, cube_length)]

        # Create particle range arrays
        if self.ndims == 3:
            x_coords = np.arange(
                particle_ranges[0][0], particle_ranges[0][1] + particle_offset_distance, particle_distance)
            y_coords = np.arange(
                particle_ranges[1][0], particle_ranges[1][1] + particle_offset_distance, particle_distance)
            z_coords = np.arange(
                particle_ranges[2][0], particle_ranges[2][1] + particle_offset_distance, particle_distance)

            # Generate the particle array using meshgrid
            xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords)
            particles = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T

        elif self.ndims == 2:
            x_coords = np.arange(
                particle_ranges[0][0], particle_ranges[0][1] + particle_offset_distance, particle_distance)
            y_coords = np.arange(
                particle_ranges[1][0], particle_ranges[1][1] + particle_offset_distance, particle_distance)

            # Generate the particle array using meshgrid for 2D
            xx, yy = np.meshgrid(x_coords, y_coords)
            particles = np.vstack((xx.ravel(), yy.ravel())).T

        else:
            raise ValueError("Only 2D and 3D cases are supported")

        # Disturb particles
        if randomness is not None:
            particles += np.random.uniform(
                -particle_offset_distance * randomness,
                particle_offset_distance * randomness,
                particles.shape)

        # Store
        self.particle_groups[self.particle_group_id] = {}
        self.particle_groups[self.particle_group_id]['particles'] = particles
        self.particle_groups[self.particle_group_id]['id'] = list(
            range(self.particles_count, self.particles_count + len(particles)))
        self.particle_groups[self.particle_group_id]['material_id'] = material_id

        # Update current particle count
        self.particles_count += len(particles)

        # Set config
        self.mpm_json["particles"].append(
            {
                "generator": {
                    "check_duplicates": True,
                    "location": f"particles_{self.particle_group_id}.txt",
                    "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                    "pset_id": self.particle_group_id,
                    "particle_type": "P3D" if self.ndims == 3 else "P2D",
                    "material_id": material_id,
                    "type": "file"
                }
            }
        )

    def add_particle_constraints(self, constraints_info:List[dict]):
        """

        Args:
            constraints_info (list[dict]): list of dicts
                Each dict contains:
                - 'pset_id': particle set id to be constrained
                - 'axis': str, one of 'x', 'y', 'z'
                - 'velocity': float, velocity

        Returns:

        """
        constraints = []

        for constraint in constraints_info:
            # Get current constraint info
            pset_id = constraint["pset_id"]
            axis = constraint["axis"]
            velocity = constraint["velocity"]

            # For each node entity set,
            for particle_set in self.entity_sets["particle_sets"]:
                # Check if current entity set is imposed a constraint by the current constraint info
                if particle_set["id"] == pset_id:
                    constraints.append({
                        "pset_id": particle_set["id"],
                        "dir": self.direction_mapping[axis],
                        "velocity": velocity
                    })
                    break

        if "boundary_conditions" not in self.mpm_json["mesh"]:
            self.mpm_json["mesh"]["boundary_conditions"] = {}

        self.mpm_json["mesh"]["boundary_conditions"]["particles_velocity_constraints"] = constraints

    def define_boundary_entity(self):
        """
        Define boundary entity set that corresponds to each boundary plain:
            entity_sets["node_sets"] = [
                {
                    "id": set_id,  # index of node set
                    "set": nodes,  # index of node
                    "axis": axis,  # x, y, z
                    "bound_loc": bound_loc  # `start` or `end`
                },
                ...,
            ]
        Note that it is hardcoded to occupy node set id from 0 to 5.
        """
        # TODO manual association of nset-id considering the
        #  `add_velocity_constraints` `add_cell_entity` `add_friction_constrains` methods.
        #   so, all the constraints' id would be better to manually assigned.
        self.entity_sets['node_sets'] = []

        # Select the appropriate axes based on the number of dimensions
        axes = AXES_3D if self.ndims == 3 else AXES_2D

        # Initialize boundary node sets
        boundary_node_ids: Dict[str, Dict[str, List[int]]] = {axis: {"start": [], "end": []} for axis in axes}

        # Populate boundary node sets
        for i, node_coord in enumerate(self.mesh_info["node_coords"]):
            for j, axis in enumerate(axes):
                if node_coord[j] == self.mesh_coord_base[j][0]:
                    boundary_node_ids[axis]["start"].append(i)
                elif node_coord[j] == self.domain_ranges[j][1]:
                    boundary_node_ids[axis]["end"].append(i)

        # Define entity set
        set_id = 0
        for axis, bounds in boundary_node_ids.items():
            for bound_loc, nodes in bounds.items():
                # TODO: may add the type of entity set that describes what this entity is about
                self.entity_sets["node_sets"].append({
                    "id": set_id,  # index of node set
                    "set": nodes,  # index of node
                    "axis": axis,  # x, y, z
                    "bound_loc": bound_loc  # `start` or `end`
                })
                set_id += 1

    def add_velocity_constraints(self, constraints_info: List[Dict]):
        """

        Args:
            constraints_info (): list of dicts
                Each dict contains:
                - 'axis': str, one of 'x', 'y', 'z'
                - 'bound_loc': str, one of 'start', 'end'
                - 'velocity': float

        Returns:

        """
        constraints = []

        for constraint in constraints_info:
            # Get current constraint info
            axis = constraint["axis"]
            bound_loc = constraint["bound_loc"]
            velocity = constraint["velocity"]

            # For each node entity set,
            for node_set in self.entity_sets["node_sets"]:
                # Check if current entity set is imposed a constraint by the current constraint info
                if node_set["axis"] == axis and node_set["bound_loc"] == bound_loc:
                    constraints.append({
                        "nset_id": node_set["id"],
                        "dir": self.direction_mapping[axis],
                        "velocity": velocity
                    })
                    break

        if "boundary_conditions" not in self.mpm_json["mesh"]:
            self.mpm_json["mesh"]["boundary_conditions"] = {}

        self.mpm_json["mesh"]["boundary_conditions"]["velocity_constraints"] = constraints

    def add_friction_constrains(self, constraints_info):
        """

        Args:
            constraints_info (): list of dicts
                Each dict contains:
                - 'axis': str, one of 'x', 'y', 'z'
                - 'bound_loc': str, one of 'start', 'end'
                - 'dir' (depreciated): int, direction index (0 for x, 1 for y, 2 for z)
                - 'sign_n': int, normal sign (-1 or 1)
                - 'friction': float

        Returns:

        """

        friction_constraints = []

        for constraint in constraints_info:
            # Get current constraint info
            axis = constraint["axis"]
            bound_loc = constraint["bound_loc"]
            sign_n = constraint["sign_n"]
            friction = constraint["friction"]

            # For each node entity set,
            for node_set in self.entity_sets["node_sets"]:
                # Check if current entity set is imposed a constraint by the current constraint info
                if node_set["axis"] == axis and node_set["bound_loc"] == bound_loc:
                    friction_constraints.append({
                        "nset_id": node_set["id"],
                        "dir": self.direction_mapping[axis],
                        "sign_n": sign_n,
                        "friction": friction
                    })
                    break

        if "boundary_conditions" not in self.mpm_json["mesh"]:
            self.mpm_json["mesh"]["boundary_conditions"] = {}

        self.mpm_json["mesh"]["boundary_conditions"]["friction_constraints"] = friction_constraints

    def add_materials(self, materials):
        """

        Args:
            materials (list): list of materials defined in dictionary.

        Returns:

        """
        self.materials = materials
        self.mpm_json["materials"] = materials

        # Check if the material model has a proper dimensionality
        for material in materials:
            if (self.ndims == 3 and "2D" in material["type"]) or (self.ndims == 2 and "3D" in material["type"]):
                raise ValueError(f"Material '{material['type']}' is not compatible with a {self.ndims}D model.")

    def add_initial_stress(
            self,
            option: str,
            top_surface: trimesh.Trimesh,
            density,
            k0=None,
            undeformed_data=None):
        """

        Args:
            top_surface (trimesh.Trimesh):
            k0 (float):
            option (str): `k0` or `stabilized_stress_data`
            undeformed_data ():

        Returns:

        """
        if self.ndims != 3:
            raise ValueError("This function currently only supports 3D domain")

        self.initial_stresses = np.zeros((self.particles_count, self.ndims))

        if option == 'k0':

            if k0 is None:
                raise ValueError("k0 is not specified")

            for pset_id, pdata in self.particle_groups.items():

                # density = find_material_property(pdata['material_id'], 'density', self.materials)

                top_zcoords = utils.get_z_coordinates(top_surface, pdata['particles'][:, [0, 1]], 'linear')
                vertical_stresses = (top_zcoords - pdata['particles'][:, 2]) * density * 9.81
                self.initial_stresses[pdata['id'], 0] = - k0 * vertical_stresses
                self.initial_stresses[pdata['id'], 1] = - k0 * vertical_stresses
                self.initial_stresses[pdata['id'], 2] = - vertical_stresses

            # raise NotImplementedError("This feature is not yet implemented")

            # if k0 is None:
            #     raise ValueError("k0 should be specified")
            #
            # gravity = self.mpm_json["external_loading_conditions"][-1]
            #
            # all_particles = []
            # for pinfo in self.mpm_json['particles']:
            #     pset_id, material_id = pinfo['generator']['pset_id'], pinfo['generator']['pset_id']
            #
            #     # Find the material associated with `material_id` and its unit weight
            #     for material_info in self.mpm_json['materials']:
            #         if material_id == material_info['id']:
            #             density = material_info['density']
            #             unit_weight = gravity * density
            #
            #     # Get the particles associated with the current `pset_id`.
            #     current_particles = self.particle_groups[pset_id]['particles']
            #
            #     # Compute stress
            #     TODO: the following won't be correct approach if the surface height is not flat but varies
            #     particle_stress = np.zeros((np.shape(current_particles)[0], 3))  # second axis is for stress xx, yy, zz
            #     particle_stress[:, 0] = k0 * (y_range[1] - particles[:, 1]) * unit_weight  # K0*H*Unit_Weight
            #     particle_stress[:, 1] = (y_range[1] - particles[:, 1]) * unit_weight  # H*Unit_Weight
            #     particle_stress[:, 2] = (y_range[1] - particles[:, 1]) * unit_weight  # H*Unit_Weight

        elif option == 'stabilized_stress_data':
            raise NotImplementedError("This feature is not yet implemented")

        else:
            raise NotImplementedError(f"Option `{option}` is invalid")

        self.mpm_json['mesh']['particles_stresses'] = 'particles_stresses.txt'

    def add_external_loadings(self, loadings):
        """

        Args:
            loadings (dict):

        Returns:

        """
        if self.ndims != len(loadings["gravity"]):
            raise ValueError("External loading does no comply with the dimension of the current simulation")

        self.mpm_json["external_loading_conditions"] = loadings

    def analysis(self, config):
        """

        Args:
            config (dict):

        Returns:

        """
        self.mpm_json["analysis"] = config

        # Check if the material model has a proper dimensionality
        if (self.ndims == 3 and "2D" in config['type']) or (self.ndims == 2 and "3D" in config['type']):
            raise ValueError(f"Analysis type '{config['type']}' is not compatible with a {self.ndims}D model.")

    def post_processing(self, config):
        """

        Args:
            config (dict):

        Returns:

        """
        self.mpm_json["post_processing"] = config

    def write(self, save_dir, file_name='mpm.json'):
        """

        Args:
            save_dir (str): directory to save the mpm input file
            file_name (str): mpm json input file name with extension

        Returns:

        """
        # Set save directory
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # --- Mesh
        print(f"Write `mesh.txt` at {save_dir}")
        mesh_file_path = os.path.join(save_dir, "mesh.txt")

        # Function to write mesh array data to file
        def append_array_to_file(file_path, array):
            with open(file_path, "a") as f:
                np.savetxt(f, array, delimiter='\t', fmt='%g')

        # Write the number of nodes and elements
        with open(mesh_file_path, "w") as f:
            f.write(f"{int(self.nnode)}\t{int(self.nele)}\n")

        # Append coordinate values of nodes and cell groups to 'mesh.txt'
        append_array_to_file(mesh_file_path, self.mesh_info["node_coords"])
        append_array_to_file(mesh_file_path, self.mesh_info["cell_groups"])

        # --- Particles
        # Function to write particle coordinates to file
        def write_coordinates_to_file(file_path, coordinates):
            with open(file_path, "w") as f:
                f.write(f"{coordinates.shape[0]} \n")
                np.savetxt(f, coordinates, delimiter='\t', fmt='%.4f')

        for pid, particle_dict in self.particle_groups.items():
            file_path = os.path.join(save_dir, f"particles_{pid}.txt")
            print(f"Write `particles_{pid}.txt` at {save_dir}")
            write_coordinates_to_file(file_path, particle_dict['particles'])

        # --- Entity
        print(f"Save `entity_sets.json`at {save_dir}")
        with open(f"{save_dir}/entity_sets.json", "w") as f:
            json.dump(self.entity_sets, f, indent=2)

        # --- particles_stresses.txt
        if self.initial_stresses is not None:
            print(f"Save `particles_stresses.txt` at {save_dir}")
            with open(f"{save_dir}/particles_stresses.txt", "w") as f:
                f.write(f"{self.initial_stresses.shape[0]}\n")
                np.savetxt(f, self.initial_stresses, delimiter='\t', fmt='%.4f')
            print('saved')

        # --- `mpm.json` config
        print(f"Save `{file_name}`at {save_dir}")
        with open(f"{save_dir}/{file_name}", "w") as f:
            json.dump(self.mpm_json, f, indent=2)

    # TODO: option to save as img
    def visualize_mesh(
            self,
            save_path,
            node_indices=False
    ):
        """
        Visualize current configuration

        Args:
            save_path (str):
            nodes (bool):
            node_indices (bool):

        Returns:

        """
        # Extract the x, y, and z coordinates
        x = self.mesh_info["node_coords"][:, 0]
        y = self.mesh_info["node_coords"][:, 1]
        z = self.mesh_info["node_coords"][:, 2]

        # Calculate the min and max values for each axis
        x_min, x_max = min(x), max(x)
        y_min, y_max = min(y), max(y)
        z_min, z_max = min(z), max(z)

        # Create a scatter plot
        fig = go.Figure(data=[go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode='markers+text' if node_indices else 'markers',
            marker=dict(
                size=3,
                opacity=0.8
            ),
            text=[str(i) for i in range(len(x))] if node_indices else None,
            textposition='top center'
        )])

        # Set plot titles and labels
        fig.update_layout(
            title='3D Mesh Plot',
            scene=dict(
                xaxis=dict(title='X Axis', range=[x_min, x_max]),
                yaxis=dict(title='Y Axis', range=[y_min, y_max]),
                zaxis=dict(title='Z Axis', range=[z_min, z_max]),
                aspectmode='data'  # Ensures equal scaling for all axes
            )
        )

        fig.write_html(save_path)
        print(f"Plot saved to {save_path}")

    def visualize_particles(
            self,
            save_path
    ):
        """

        Args:
            save_path (str):

        Returns:

        """
        if self.ndims == 2:
            raise ValueError("This feature does not support 2D case yet.")

        # Create a figure
        fig = go.Figure()

        # Define color palette
        # colors = ['Viridis', 'Cividis', 'Plasma', 'Inferno', 'Magma', 'Turbo']

        for i, particle_dict in self.particle_groups.items():
            # Extract the x, y, and z coordinates
            x = particle_dict['particles'][:, 0]
            y = particle_dict['particles'][:, 1]
            z = particle_dict['particles'][:, 2]

            # Create a scatter plot for each group
            fig.add_trace(go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(
                    size=3,
                    # color=colors[i],
                    opacity=0.8
                ),
                name=f'Group-{i}'
            ))

        # Set plot titles and labels
        fig.update_layout(
            title='3D Particles Plot',
            scene=dict(
                xaxis=dict(title='X Axis', range=self.domain_ranges[0]),
                yaxis=dict(title='Y Axis', range=self.domain_ranges[1]),
                zaxis=dict(title='Z Axis', range=self.domain_ranges[2]),
                aspectmode='data'  # Ensures equal scaling for all axes
            )
        )

        fig.write_html(save_path)
        print(f"Plot saved to {save_path}")


def find_material_property(id, field, material_list):
    for item in material_list:
        if item['id'] == id:
            return item[field]


def get_h5(directory, timestep, mpi):
    # Create an empty list to store DataFrames
    dfs = []

    # Iterate over different files from MPI and append to list
    for i in range(mpi):
        file = f'particles-{i}_{mpi}-{timestep}.h5'  # ex) particles-26_32-0120000.h5
        h5_path = os.path.join(directory, file)
        dfs.append(pd.read_hdf(h5_path, 'table'))

    # Concatenate all DataFrames
    df = pd.concat(dfs, ignore_index=True)

    return df
