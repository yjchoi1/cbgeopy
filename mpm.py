import numpy as np
import plotly.graph_objects as go
import pandas as pd
import os
import json
import trimesh
import utils
import random
import argparse


class MPMConfig:
    def __init__(
            self,
            domain_origin,
            domain_length,
            title="mpm_input_config"
    ):
        # Mesh
        self.cell_size = None
        self.nnode = None
        self.nele = None
        self.mesh_info = None
        self.mesh_coord_base = None

        # Entity
        self.entity_sets = None

        # Particle
        self.particle_group_id = 0
        self.particle_groups = {}

        # General
        self.direction_mapping = {"x": 0, "y": 1, "z": 2}
        self.ndims = len(domain_length)
        self.domain_origin = domain_origin
        self.domain_length = domain_length
        if not all(len(lst) == 3 for lst in [domain_origin, domain_length]):
            raise NotImplementedError("Current version only support 3D")
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

    def add_mesh(self, cell_size):
        """
        Make mesh coordinate array & cell group
        Args:
            cell_size (list): [x_size, y_size, z_size]

        Returns:

        """
        self.cell_size = cell_size

        # Error handling
        if not all(i == cell_size[0] for i in cell_size):
            raise NotImplementedError("Current version only support the same element size for all dimension")

        # Make coordinate base for each axis
        mesh_coord_base = []
        for dim, domain_range in enumerate(self.domain_ranges):
            coord_base = np.arange(
                domain_range[0], domain_range[1] + cell_size[dim], cell_size[dim])
            mesh_coord_base.append(coord_base)
        self.mesh_coord_base = mesh_coord_base

        # Create node coordinates
        node_coords = []
        for z in self.mesh_coord_base[2]:
            for y in self.mesh_coord_base[1]:
                for x in self.mesh_coord_base[0]:
                    node_coords.append([x, y, z])
        node_coords = np.array(node_coords)

        # Compute the number of nodes and elements (cells) for each dimension
        # - node
        nnode_x = len(mesh_coord_base[0])
        nnode_y = len(mesh_coord_base[1])
        nnode_z = len(mesh_coord_base[2])
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
            particle_group_id (int):
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
        self.particle_groups[self.particle_group_id] = particles

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

    def add_particles_from_topology(
            self,
            lower_topology: trimesh.Trimesh,
            upper_topology: trimesh.Trimesh,
            n_particle_per_cell: int,
            material_id: int,
            particle_group_id: int = None
    ):

        # Assign a particle group id
        if particle_group_id is None:
            self.particle_group_id += 1
        else:
            self.particle_group_id = particle_group_id

        # Fill particles between two meshes
        particles = utils.fill_particles_between_mesh(
            lower_topology, upper_topology, self.cell_size, n_particle_per_cell)

        # Store
        self.particle_groups[self.particle_group_id] = particles

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

    def add_particles_cube(
            self,
            cube_origin,
            cube_length,
            material_id,
            n_particle_per_cell,
            particle_group_id=None):
        """

        Args:
            cube_origin (list): the lower edge of the cube, [x_min, y_min, z_min]
            cube_length (list): the length of the cube defined from `cube_origin` [x_len, y_len, z_len]
            material_id (int):
            n_particle_per_cell (int):
            particle_group_id (int):

        Returns:

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
        x_coords = np.arange(
            particle_ranges[0][0], particle_ranges[0][1] + particle_offset_distance, particle_distance)
        y_coords = np.arange(
            particle_ranges[1][0], particle_ranges[1][1] + particle_offset_distance, particle_distance)
        z_coords = np.arange(
            particle_ranges[2][0], particle_ranges[2][1] + particle_offset_distance, particle_distance)

        # Generate the particle array using meshgrid
        xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords)
        particles = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T

        # Store
        self.particle_groups[self.particle_group_id] = particles

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

    # TODO: k0 stress

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
        """
        self.entity_sets = {"node_sets": []}

        # Define boundary node sets
        boundary_node_ids = {axis: {"start": [], "end": []} for axis in ["x", "y", "z"]}

        # Get boundaries
        x_bounds, y_bounds, z_bounds = self.domain_ranges

        # Node boundary entity
        for i, node_coord in enumerate(self.mesh_info["node_coords"]):
            if node_coord[0] == x_bounds[0]:
                boundary_node_ids["x"]["start"].append(i)
            if node_coord[0] == x_bounds[1]:
                boundary_node_ids["x"]["end"].append(i)
            if node_coord[1] == y_bounds[0]:
                boundary_node_ids["y"]["start"].append(i)
            if node_coord[1] == y_bounds[1]:
                boundary_node_ids["y"]["end"].append(i)
            if node_coord[2] == z_bounds[0]:
                boundary_node_ids["z"]["start"].append(i)
            if node_coord[2] == z_bounds[1]:
                boundary_node_ids["z"]["end"].append(i)

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

    def add_velocity_constraints(self, constraints_info):
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
                - 'dir': int, direction index (0 for x, 1 for y, 2 for z)
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
        self.mpm_json["materials"] = materials

    def add_external_loadings(self, loadings):
        """

        Args:
            loadings (dict):

        Returns:

        """
        self.mpm_json["external_loading_conditions"] = loadings

    def analysis(self, config):
        """

        Args:
            config (dict):

        Returns:

        """
        self.mpm_json["analysis"] = config

    def post_processing(self, config):
        """

        Args:
            config (dict):

        Returns:

        """
        self.mpm_json["post_processing"] = config

    def write(self, save_dir):
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

        for pid, particle_coords in self.particle_groups.items():
            file_path = os.path.join(save_dir, f"particles_{pid}.txt")
            print(f"Write `particles_{pid}.txt` at {save_dir}")
            write_coordinates_to_file(file_path, particle_coords)

        # --- Entity
        print(f"Save `entity_sets.json`at {save_dir}")
        with open(f"{save_dir}/entity_sets.json", "w") as f:
            json.dump(self.entity_sets, f, indent=2)

        # --- `mpm.json` config
        print(f"Save `mpm.json`at {save_dir}")
        with open(f"{save_dir}/mpm.json", "w") as f:
            json.dump(self.mpm_json, f, indent=2)

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
                aspectmode='cube'  # Ensures equal scaling for all axes
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

        # Create a figure
        fig = go.Figure()

        # Define color palette
        # colors = ['Viridis', 'Cividis', 'Plasma', 'Inferno', 'Magma', 'Turbo']

        for i, particles in self.particle_groups.items():
            # Extract the x, y, and z coordinates
            x = particles[:, 0]
            y = particles[:, 1]
            z = particles[:, 2]

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







