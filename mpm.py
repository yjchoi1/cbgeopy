"""MPM configuration module for setting up Material Point Method simulations.

This module provides the MPMConfig class and helper functions for configuring MPM simulations.
It handles mesh generation, particle placement, material properties, boundary conditions and more.
"""

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
import gstools as gs
from scipy.spatial import cKDTree, KDTree
from scipy import interpolate
from typing import List, Dict, Callable, Optional
import shapely
import shapely.geometry as geom
from cbgeopy.utils import generate_random_field


# Define constants for axes
AXES_3D = ["x", "y", "z"]
AXES_2D = ["x", "y"]


class MPMConfig:
    """Class for configuring Material Point Method simulations.

    This class handles all aspects of MPM simulation setup including:
    - Mesh generation and configuration
    - Particle placement and properties 
    - Material definitions
    - Boundary conditions
    - Analysis settings
    - Output configuration

    Attributes:
        initial_stresses: Initial stress states for particles
        materials: List of material definitions
        cell_size: Size of mesh cells
        n_cells_per_dim: Number of cells per dimension
        nnode: Total number of nodes
        nele: Total number of elements
        mesh_info: Dictionary containing mesh information
        mesh_coord_base: Base coordinates for mesh generation
        entity_sets: Dictionary of entity sets
        particle_group_id: Current particle group ID
        particles_count: Total number of particles
        particle_groups: Dictionary of particle groups
        direction_mapping: Maps direction names to indices
        ndims: Number of dimensions
        domain_origin: Origin coordinates of domain
        domain_length: Length of domain in each dimension
        domain_ranges: Ranges of domain in each dimension
        mpm_json: Dictionary containing MPM configuration
    """

    def __init__(
            self,
            domain_origin: List[float],
            domain_length: List[float],
            title: str = "mpm_input_config"
    ):
        """Initialize MPM configuration.

        Args:
            domain_origin: Origin coordinates of simulation domain
            domain_length: Length of domain in each dimension
            title: Title for the configuration
        """
        # Mesh
        self.initial_stresses = None
        self.materials = []
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
        """Make mesh coordinate array & cell group.

        Args:
            n_cells_per_dim: Number of cells per dimension [nx, ny, nz]
            outer_cell_thickness: Thickness of outer mesh layer. If provided, adds outer mesh
                beyond domain_ranges. If 0, mesh matches domain_ranges exactly.

        Raises:
            ValueError: If outer_cell_thickness is negative
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
        """Add particles from a CSV file.

        Args:
            path: Path to CSV file containing particle coordinates
            material_id: Material ID to associate with these particles
            particle_group_id: Optional particle group ID, auto-increments if None

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
        
    def add_particles_from_random_field(
        self,
        polygons_params: List,
        n_particle_per_cell: int,
        randomness: float = None,
        random_field_method: str = "gstools"
    ):
        """Add particles within polygons with random field values.

        Args:
            polygons_params: List of dictionaries containing:
                - polygon_points: List of points defining polygon vertices
                - random_params: Dict with random field parameters (mean, std, len_scale)
                - custom_correlation: (Optional) Custom correlation function (gstools.CovModel)
            n_particle_per_cell: Number of particles per cell per dimension
            randomness: Factor for random perturbation of particle positions
            random_field_method: Method to generate random field ("gstools" or "utils")

        Note:
            Unlike other add_particles methods, this automatically sets particle_group_id based
            on the particles' associated cell id, and material_id equals particle_group_id.

        Raises:
            ValueError: If used with 3D domain or if no particles found in polygon
        """
        if self.ndims == 3:
            raise ValueError("This feature is only for 2D domain")

        # Generate particle grid
        candidate_particles, _, particle_offset_distance = self._generate_particle_grid(
            self.domain_origin,
            self.domain_length,
            self.cell_size,
            n_particle_per_cell,
            self.ndims
        )
        
        # Get mesh grid coordinates
        x = self.mesh_coord_base[0]  # node coordinates in x direction
        y = self.mesh_coord_base[1]  # node coordinates in y direction
        
        # make cell-particle mapping
        self.cell_particle_field_groups = []
        
        for region in polygons_params:
            # Validate required keys
            if "polygon_points" not in region:
                raise ValueError("Missing required 'polygon_points' in polygon parameters")
            if "random_params" not in region:
                raise ValueError("Missing required 'random_params' in polygon parameters")
                
            # Validate random_params structure
            required_params = ["mean", "std", "len_scale"]
            for param in required_params:
                if param not in region["random_params"]:
                    raise ValueError(f"Missing required '{param}' in random_params")
            
            poly = shapely.geometry.Polygon(region["polygon_points"])
            random_params = region["random_params"]
            
            # Generate random field for cells
            if "custom_correlation" in region and region["custom_correlation"] is not None:
                # Create instance of custom correlation and set parameters
                model = region["custom_correlation"]
                model.dim = self.ndims
                model.var = random_params["std"]**2
                model.len_scale = random_params["len_scale"]
                
                # Use GSTools for custom correlation
                srf = gs.SRF(model, mean=random_params["mean"])
                # Generate values for cell corners
                field_values = srf.structured((x[:-1], y[:-1]))
            else:
                # Choose between GSTools and utils.generate_random_field
                if random_field_method.lower() == "gstools":
                    # Use GSTools with Gaussian model
                    model = gs.Gaussian(
                        dim=self.ndims, 
                        var=random_params["std"]**2, 
                        len_scale=random_params["len_scale"]
                    )
                    srf = gs.SRF(model, mean=random_params["mean"])
                    # Generate values for cell corners
                    field_values = srf.structured((x[:-1], y[:-1]))
                elif random_field_method.lower() == "utils":                    
                    # Generate the random field using the function from utils
                    field_values, _, _ = generate_random_field(
                        x_range=(x[0], x[-2]),
                        y_range=(y[0], y[-2]),
                        nx=len(x) - 1,
                        ny=len(y) - 1,
                        lx=random_params["len_scale"][0] if isinstance(random_params["len_scale"], (list, tuple)) else random_params["len_scale"],
                        ly=random_params["len_scale"][1] if isinstance(random_params["len_scale"], (list, tuple)) else random_params["len_scale"],
                        sigma2=random_params["std"]**2,
                        mean_value=random_params["mean"],
                        seed=random_params.get("seed", None)
                    )
                else:
                    raise ValueError(f"Invalid random_field_method: {random_field_method}. Choose 'gstools' or 'utils'")
            
            # Create points for all candidate particles
            points = [shapely.geometry.Point(p) for p in candidate_particles]
            
            # Filter particles inside polygon
            mask = [poly.contains(point) for point in points]
            particles = candidate_particles[mask]
            
            if len(particles) == 0:
                raise ValueError("No particles found within the polygon")
                                
            # Disturb particles if randomness is specified
            if randomness is not None:
                particles += np.random.uniform(
                    -particle_offset_distance * randomness,
                    particle_offset_distance * randomness,
                    particles.shape
                )
                
            ################## Start to make cell-particle-field mapping ##################
            # We need to map which particles are associated with which cell since 
            #   the randomfield values are cell-based. If the particles in a certain cell, the corresponding
            #   material property should follow the cell's field value. 
            
            # Find cell indices for each particle
            particle_cell_ids_x = np.searchsorted(x, particles[:, 0]) - 1
            particle_cell_ids_y = np.searchsorted(y, particles[:, 1]) - 1
            
            # Get field values for each particle's cell
            particle_field_values = field_values[particle_cell_ids_x, particle_cell_ids_y]
            
            # Convert to linear cell indices
            particle_cell_ids = particle_cell_ids_y * (len(x) - 1) + particle_cell_ids_x
            
            unique_cells = np.unique(particle_cell_ids)
            
            # Define cell-particle-field mapping and make materials for each cell
            for cell_id in unique_cells:
                # Get indices of particles in this cell
                particle_mask = (particle_cell_ids == cell_id)
                particle_ids = np.where(particle_mask)[0]
                current_particles = np.array(particles[particle_mask])
                
                # Get field values for particles in this cell
                cell_field_values = particle_field_values[particle_mask]
                
                # Store the cell-particle mapping
                self.cell_particle_field_groups.append({
                    "cell_id": int(cell_id),
                    "particle_ids": particle_ids.tolist(),
                    "particles": current_particles,
                    "field_values": cell_field_values.tolist()
                })
                
            ################## End of cell-particle-field mapping ##################

        # Set config
        for cpf_map in self.cell_particle_field_groups:
            # Before assigning particle groups, find the next available ID
            if len(self.particle_groups) == 0:
                next_id = 0
            else:
                next_id = max(self.particle_groups.keys()) + 1

            # Use this ID directly without adding cell_id
            self.particle_group_id = next_id
                
            # Store particles and field values
            self.particle_groups[self.particle_group_id] = {}
            self.particle_groups[self.particle_group_id]['particles'] = cpf_map["particles"]
            self.particle_groups[self.particle_group_id]['id'] = cpf_map["particle_ids"]
            self.particle_groups[self.particle_group_id]['field_values'] = cpf_map["field_values"]
            
            # Update current particle count
            self.particles_count += len(cpf_map["particles"])
            
            # Set particle generator config
            self.mpm_json["particles"].append(
                {
                    "generator": {
                        "check_duplicates": True,
                        "location": f"particles_{self.particle_group_id}.txt",
                        "io_type": "Ascii3D" if self.ndims == 3 else "Ascii2D",
                        "pset_id": self.particle_group_id,
                        "particle_type": "P3D" if self.ndims == 3 else "P2D",
                        "material_id": self.particle_group_id,
                        "type": "file"
                    }
                }
            )
            
    def add_particles_from_polygon(
        self,
        polygon_info: List,
        n_particle_per_cell: int,
        randomness: float = None
    ):
        """Add particles within defined polygons.

        Args:
            polygon_info: List of dictionaries containing:
                - polygon_points: List of points defining polygon vertices
                - material_id: Material ID for particles in this polygon
                - particle_group_id: (Optional) Group ID for these particles
            n_particle_per_cell: Number of particles per cell per dimension
            randomness: Factor for random perturbation of particle positions

        Raises:
            ValueError: If used with 3D domain
        """
        if self.ndims == 3:
            raise ValueError("This feature is only for 2D domain")

        # Generate particle grid
        candidate_particles, _, particle_offset_distance = self._generate_particle_grid(
            self.domain_origin,
            self.domain_length,
            self.cell_size,
            n_particle_per_cell,
            self.ndims
        )

        # Process each polygon
        for poly in polygon_info:
            # Assign a particle group id
            if "particle_group_id" in poly:
                self.particle_group_id = poly["particle_group_id"]
            else:
                self.particle_group_id += 1

            # Create Shapely polygon
            polygon = shapely.geometry.Polygon(poly["polygon_points"])
            
            # Create Shapely points for all candidate particles
            points = [shapely.geometry.Point(p) for p in candidate_particles]
            
            # Filter particles inside polygon
            mask = [polygon.contains(point) for point in points]
            particles = candidate_particles[mask]

            # Disturb particles if randomness is specified
            if randomness is not None and len(particles) > 0:
                particles += np.random.uniform(
                    -particle_offset_distance * randomness,
                    particle_offset_distance * randomness,
                    particles.shape)

            # Store particles
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
                        "io_type": "Ascii2D",
                        "pset_id": self.particle_group_id,
                        "particle_type": "P2D",
                        "material_id": poly["material_id"],
                        "type": "file"
                    }
                }
            )

    def add_particles_from_lines(
            self,
            layer_info: List,
            n_particle_per_cell: int,
    ):
        """Add particles in layers defined by line segments.

        Args:
            layer_info: List of dictionaries containing:
                - line_points: Points defining upper boundary line of layer
                - material_id: Material ID for this layer
                - particle_group_id: (Optional) Group ID for these particles
                - randomness: (Optional) Factor for particle position randomization
            n_particle_per_cell: Number of particles per cell per dimension

        Raises:
            ValueError: If used with 3D domain

        Note:
            Line points should span entire x-domain range.
            First layer assumes y=0 as lower boundary.
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
            if "randomness" in layer:
                particles += np.random.uniform(
                    -particle_offset_distance * layer["randomness"],
                    particle_offset_distance * layer["randomness"],
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
        """Add particles between two topographic surfaces.

        This method fills the space between two topographic surfaces with particles. The particles
        are placed in a regular grid pattern and can optionally be randomly perturbed.

        Args:
            lower_topography: Mesh defining the lower surface topography
            upper_topography: Mesh defining the upper surface topography 
            n_particle_per_cell: Number of particles per cell per dimension
            material_id: ID of material to assign to this particle set
            z_find_method: Method to find z-coordinate of mesh ('' or '')
            base_find_method: Method to find base area where surfaces overlap
                'alphashape': Estimates x-y projection area using alphashape
                'simple': Uses max/min x,y values, restricted to square shape
            z_fill_method: Method to fill between lower and upper z-coordinates
                'simple': Uses exact z-coordinate of mesh
                'round': Uses nearest particle grid points
            randomness: Factor to randomly perturb particle positions
            overlap_tolerance: Maximum distance to consider particles as overlapping
            particle_group_id: ID to assign to this particle group

        Raises:
            ValueError: If used with 2D domain (only supports 3D)

        Note:
            The particles are stored in self.particle_groups with the specified particle_group_id.
            If particle_group_id is None, it will be auto-incremented.
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
        """Define particle entity sets for the MPM solver.

        This method creates particle entity sets that define different materials for each particle set.
        The entity sets are used as input to the MPM solver to specify material properties.
        Each particle set is assigned a unique ID and contains a list of particle indices.

        Example:
            For two particle groups with indices [0,1,2] and [3,4,5], the entity sets would be::

                {
                    'particle_sets': [
                        {'id': 0, 'set': [0,1,2]},
                        {'id': 1, 'set': [3,4,5]} 
                    ]
                }

        Note:
            The particle sets are stored in self.entity_sets['particle_sets'] as a list of
            dictionaries. Each dictionary contains:
            
            - id: Unique identifier for the particle set (int)
            - set: List of particle indices belonging to this set (List[int])

            For more details on entity set format, see:
            https://mpm.cb-geo.com/#/user/preprocess/entity-sets
        """
        self.entity_sets['particle_sets'] = []
        for set_id, particle_dict in self.particle_groups.items():
            self.entity_sets["particle_sets"].append({
                "id": set_id,  # index of particle set
                "set": particle_dict['id'],  # index of particles for current set
            })
                
    def remove_overlapping_particles(self, overlap_tolerance: float):
        """Remove overlapping particles between different particle groups.
        
        This method iterates through particle groups in order and removes any particles 
        that overlap with particles from previous groups. A particle is considered overlapping
        if it is within the overlap_tolerance distance of any existing particle.
        
        The particle IDs are reordered after removing overlaps to maintain consecutive numbering.
        
        Args:
            overlap_tolerance: Maximum distance between particles to consider them as overlapping.
                             Particles closer than this distance will be considered overlapping
                             and one will be removed.
        
        Raises:
            ValueError: If no particle groups exist
            
        Note:
            This modifies the particle_groups dictionary in place, updating both the particles
            and their IDs for each group.
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
        """Add a cell entity set based on coordinate ranges.

        Args:
            nset_id: Node set ID to assign to this entity set. Must be > 5 since 0-5 are reserved.
            ranges: List of coordinate ranges [[x_min, x_max], [y_min, y_max], [z_min, z_max]]

        Raises:
            ValueError: If nset_id <= 5 or if ranges length doesn't match domain dimensions

        Note:
            The cell entity set is stored in self.entity_sets["node_sets"] and contains all nodes
            within the specified coordinate ranges.
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
        """Add a group of particles within a defined cubic region.

        Args:
            cube_origin: Lower edge coordinates of cube [x_min, y_min, (z_min)]
            cube_length: Lengths of cube along each axis [x_len, y_len, (z_len)]
            material_id: Material ID to assign to this particle group
            n_particle_per_cell: Number of particles per dimension in each cell
            randomness: Factor to randomly perturb particle positions
            particle_group_id: ID to assign to this particle group

        Raises:
            ValueError: If domain dimensions are not 2D or 3D

        Note:
            The particles are stored in self.particle_groups with the specified particle_group_id.
            If particle_group_id is None, it will be auto-incremented.
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
            range(self.particles_count, self.particles_count + len(particles))
        )
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
            constraints_info: List of constraint dictionaries, each containing:
                - 'pset_id': Particle set ID to constrain
                - 'axis': Direction to constrain ('x', 'y', or 'z')
                - 'velocity': Velocity value to apply

        Note:
            The constraints are added to self.mpm_json["mesh"]["boundary_conditions"]
            under the "particles_velocity_constraints" key.
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
        """Define boundary entity sets for each boundary plane.

        Creates node sets for each boundary plane of the domain. The node sets are stored in
        self.entity_sets["node_sets"] with IDs 0-5 reserved for boundaries.

        Each node set dictionary contains:
            - id: Node set ID (0-5)
            - set: List of node indices
            - axis: Boundary axis ('x', 'y', 'z')
            - bound_loc: Boundary location ('start' or 'end')

        Note:
            Node set IDs 0-5 are hardcoded for boundaries. Other entity sets should use IDs > 5.
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
        """Add velocity constraints to boundary nodes.

        Args:
            constraints_info: List of constraint dictionaries, each containing:
                - 'axis': Direction to constrain ('x', 'y', or 'z')
                - 'bound_loc': Boundary location ('start' or 'end')
                - 'velocity': Velocity value to apply

        Note:
            The constraints are added to self.mpm_json["mesh"]["boundary_conditions"]
            under the "velocity_constraints" key.
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
        """Add friction constraints to boundary nodes.

        Args:
            constraints_info: List of constraint dictionaries, each containing:
                - 'axis': Direction to constrain ('x', 'y', or 'z')
                - 'bound_loc': Boundary location ('start' or 'end')
                - 'sign_n': Normal direction sign (-1 or 1)
                - 'friction': Friction coefficient

        Note:
            The constraints are added to self.mpm_json["mesh"]["boundary_conditions"]
            under the "friction_constraints" key.
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

    def add_materials(
        self, 
        materials: Optional[List[Dict]] = None, *, 
        option: Optional[str] = None,
        material_type: Optional[str] = None
    ):
        """Add materials to the MPM simulation.
        
        Args:
            materials: List of material dictionaries with properties. Required unless using a special option.
            option: Special material generation option. One of: {None, "random_field"}
            material_type: Type of material model. One of: {None, "MohrCoulomb2D"}

        Special Options:
            random_field: Generates materials with random field values for each cell particle group.
                Requires material_type="MohrCoulomb2D"  
            Random field generation:
                mpm.add_materials(option="random_field", material_type="MohrCoulomb2D")

        Raises:
            ValueError: If arguments are invalid or incompatible with model dimensions
        """
        # Define valid options
        VALID_OPTIONS = {None, "random_field"}
        VALID_TYPES = {None, "MohrCoulomb2D"}
        
        # Validate inputs
        if option not in VALID_OPTIONS:
            raise ValueError(f"Invalid option: {option}. Must be one of: {VALID_OPTIONS}")
        if material_type not in VALID_TYPES:
            raise ValueError(f"Invalid material_type: {material_type}. Must be one of: {VALID_TYPES}")
        
        # Handle special options
        if option == "random_field":
            if material_type != "MohrCoulomb2D":
                raise ValueError(
                    "random_field option requires material_type='MohrCoulomb2D'")
            if not hasattr(self, 'particle_groups'):
                raise ValueError(
                    "random_field option requires particle_groups to be defined")
            
            # Generate materials for each cell group
            materials = []
            for group_id, particle_group_info in self.particle_groups.items():
                if "field_values" in particle_group_info:
                    materials.append({
                        "id": group_id,
                        "type": material_type,
                        "density": 1800,
                        "youngs_modulus": 40000000.0,
                        "poisson_ratio": 0.3,
                        "friction": particle_group_info["field_values"][0],
                        "dilation": 0.0,
                        "cohesion": 1000,
                        "tension_cutoff": 100,
                        "softening": False,
                        "peak_pdstrain": 0.0,
                        "residual_friction": 30.0,
                        "residual_dilation": 0.0,
                        "residual_cohesion": 0.0,
                        "residual_pdstrain": 0.0
                    })
        
        # Standard material definition
        else:
            if not materials:
                raise ValueError(
                    "materials argument is required when not using a special option")
        
        # Validate material dimensionality
        for material in materials:
            if "type" not in material:
                raise ValueError(f"Material missing required 'type' field: {material}")
            if (self.ndims == 3 and "2D" in material["type"]) or (self.ndims == 2 and "3D" in material["type"]):
                raise ValueError(f"Material '{material['type']}' is not compatible with a {self.ndims}D model.")
        
        # Store the materials
        if not hasattr(self, 'materials') or self.materials is None:
            self.materials = []
        self.materials.extend(materials)
        
        # Update MPM JSON config
        if "materials" not in self.mpm_json:
            self.mpm_json["materials"] = []
        self.mpm_json["materials"].extend(materials)

    def add_initial_stress(
            self,
            option: str,
            top_surface: trimesh.Trimesh,
            density,
            k0=None,
            undeformed_data=None):
        """Add initial stress state to particles.

        Args:
            option: Method to calculate initial stress ('k0' or 'stabilized_stress_data')
            top_surface: Mesh defining the top surface topography
            density: Material density
            k0: Lateral earth pressure coefficient (required if option='k0')
            undeformed_data: Data for stabilized stress calculation (required if option='stabilized_stress_data')

        Raises:
            ValueError: If used with 2D domain or if required parameters are missing
            NotImplementedError: If option is not implemented

        Note:
            The initial stresses are stored in self.initial_stresses and written to
            'particles_stresses.txt' when the configuration is saved.
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
        """Add external loading conditions.

        Args:
            loadings: Dictionary containing loading conditions, including:
                - gravity: List of gravity components matching domain dimensions

        Raises:
            ValueError: If loading dimensions don't match domain dimensions

        Note:
            The loadings are stored in self.mpm_json["external_loading_conditions"].
        """
        if self.ndims != len(loadings["gravity"]):
            raise ValueError("External loading does no comply with the dimension of the current simulation")

        self.mpm_json["external_loading_conditions"] = loadings

    def analysis(self, config):
        """Configure analysis settings for the MPM simulation.

        Args:
            config (dict): Dictionary containing analysis configuration parameters including:
                - type: Analysis type (e.g. 'MPMExplicit2D', 'MPMExplicit3D')
                - Other analysis-specific parameters

        Raises:
            ValueError: If analysis type dimensionality doesn't match domain dimensions

        Note:
            The analysis settings are stored in self.mpm_json["analysis"].
            Analysis type must match domain dimensionality (2D/3D).
        """
        self.mpm_json["analysis"] = config

        # Check if the material model has a proper dimensionality
        if (self.ndims == 3 and "2D" in config['type']) or (self.ndims == 2 and "3D" in config['type']):
            raise ValueError(f"Analysis type '{config['type']}' is not compatible with a {self.ndims}D model.")

    def post_processing(self, config):
        """Configure post-processing settings for the MPM simulation.

        Args:
            config (dict): Dictionary containing post-processing configuration parameters.

        Note:
            The post-processing settings are stored in self.mpm_json["post_processing"].
        """
        self.mpm_json["post_processing"] = config

    def write(self, save_dir, file_name='mpm.json', save_options=None):
        """Write MPM configuration and data files to disk.

        Args:
            save_dir (str): Directory path to save the MPM input files
            file_name (str): Name of the MPM JSON input file with extension
            save_options (dict, optional): Dictionary controlling which files to save.
                Keys are file types, values are booleans. Available options:
                {
                    'mesh': True,           # mesh.txt
                    'particles': True,      # particles_*.txt
                    'entity_sets': True,    # entity_sets.json
                    'stresses': True,       # particles_stresses.txt
                    'config': True          # mpm.json
                }
                If None, all files will be saved.

        Note:
            This method writes several files based on save_options:
            - mesh.txt: Contains mesh node coordinates and cell groups
            - particles_*.txt: Particle coordinates for each particle group
            - entity_sets.json: Entity set definitions
            - particles_stresses.txt: Initial particle stresses if defined
            - mpm.json: Main configuration file
        """
        # Default save options
        default_options = {
            'mesh': True,
            'particles': True,
            'entity_sets': True,
            'stresses': True,
            'config': True
        }
        
        # Use provided options or defaults
        save_options = save_options or default_options
        save_options = {**default_options, **save_options}  # Merge with defaults

        # Set save directory
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # --- Mesh
        if save_options['mesh']:
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
        if save_options['particles']:
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
        if save_options['entity_sets']:
            print(f"Save `entity_sets.json`at {save_dir}")
            with open(f"{save_dir}/entity_sets.json", "w") as f:
                json.dump(self.entity_sets, f, indent=2)

        # --- particles_stresses.txt
        if save_options['stresses'] and self.initial_stresses is not None:
            print(f"Save `particles_stresses.txt` at {save_dir}")
            with open(f"{save_dir}/particles_stresses.txt", "w") as f:
                f.write(f"{self.initial_stresses.shape[0]}\n")
                np.savetxt(f, self.initial_stresses, delimiter='\t', fmt='%.4f')
            print('saved')

        # --- `mpm.json` config
        if save_options['config']:
            print(f"Save `{file_name}`at {save_dir}")
            with open(f"{save_dir}/{file_name}", "w") as f:
                json.dump(self.mpm_json, f, indent=2)

    # TODO: option to save as img
    def visualize_mesh(
            self,
            save_path,
            node_indices=False
    ):
        """Visualize the mesh configuration in 3D.

        Args:
            save_path (str): Path to save the visualization HTML file
            node_indices (bool): Whether to display node indices in the visualization

        Note:
            Creates an interactive 3D plot using Plotly showing mesh nodes.
            The plot is saved as an HTML file at the specified path.
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
        """Visualize particle groups in 3D.

        Args:
            save_path (str): Path to save the visualization HTML file

        Raises:
            ValueError: If attempting to visualize 2D particle configuration

        Note:
            Creates an interactive 3D plot using Plotly showing particles from all groups.
            Each particle group is plotted in a different color.
            The plot is saved as an HTML file at the specified path.
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

    @staticmethod
    def _generate_particle_grid(
        domain_origin, 
        domain_length, 
        cell_size, 
        n_particle_per_cell, 
        ndims
        ):
        """Generate a uniform grid of particle coordinates.
        
        Args:
            domain_origin (List[float]): Origin coordinates of the domain
            domain_length (List[float]): Length of domain in each dimension
            cell_size (List[float]): Size of cells in each dimension
            n_particle_per_cell (int): Number of particles per cell per dimension
            ndims (int): Number of dimensions (2 or 3)
            
        Returns:
            tuple: Contains:
                - particles (np.ndarray): Array of particle coordinates
                - particle_distance (float): Distance between particles
                - particle_offset_distance (float): Offset distance from domain boundaries
        """
        particle_distance = cell_size[0] / n_particle_per_cell
        particle_offset_distance = particle_distance / 2
        
        # Create particle range arrays that cover the whole domain
        particle_ranges = [
            (origin + particle_offset_distance, origin + length - particle_offset_distance)
            for origin, length in zip(domain_origin, domain_length)
        ]
        
        # Generate coordinate arrays for each dimension
        coords = [
            np.arange(prange[0], prange[1] + particle_offset_distance, particle_distance)
            for prange in particle_ranges[:ndims]
        ]
        
        # Generate particle grid using meshgrid
        if ndims == 2:
            xx, yy = np.meshgrid(coords[0], coords[1])
            particles = np.vstack((xx.ravel(), yy.ravel())).T
        else:  # ndims == 3
            xx, yy, zz = np.meshgrid(coords[0], coords[1], coords[2])
            particles = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
        
        return particles, particle_distance, particle_offset_distance


def find_material_property(id, field, material_list):
    """Find a specific property of a material by its ID.

    Args:
        id (int): Material ID to search for
        field (str): Name of the material property to retrieve
        material_list (List[dict]): List of material dictionaries

    Returns:
        The value of the specified field for the material with matching ID
    """
    for item in material_list:
        if item['id'] == id:
            return item[field]


def get_h5(directory, timestep, mpi):
    """Read and combine HDF5 particle data files from MPI processes.

    Args:
        directory (str): Directory containing the HDF5 files
        timestep (int): Timestep to read
        mpi (int): Number of MPI processes

    Returns:
        pandas.DataFrame: Combined particle data from all MPI processes
    """
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
