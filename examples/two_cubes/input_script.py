import trimesh
import utils
from mpm import MPMConfig
import demo_utils
from functools import partial
import numpy as np
import vis_utils
import os
import trimesh


save_dir = 'examples/two_cubes/'

# Set config
lx, ly, lz = 1.0, 1.0, 1.0
mpm = MPMConfig(domain_origin=[0, 0, 0], domain_length=[lx, ly, lz])

# Mesh
cell_size = 0.2
mpm.add_mesh(
    n_cells_per_dim=[int(lx/cell_size), int(ly/cell_size), int(lz/cell_size)])

# Add materials
mpm.add_materials(
    [
        {
            "id": 0,
            "density": 1800,
            "youngs_modulus": 2000000.0,
            "poisson_ratio": 0.3,
            "friction": 30.0,
            "dilation": 0.0,
            "cohesion": 100,
            "tension_cutoff": 50,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 10.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb3D"
        },
        {
            "id": 1,
            "density": 1800,
            "youngs_modulus": 20000000.0,
            "poisson_ratio": 0.35,
            "friction": 45.0,
            "dilation": 0.0,
            "cohesion": 20e6,
            "tension_cutoff": 50,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 45.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb3D"
        }
    ]
)

# Particles
mpm.add_particles_cube(
    cube_origin=[0, 0, 0],
    cube_length=[0.4, 0.4, 0.4],
    material_id=0,
    n_particle_per_cell=4,
    particle_group_id=0
)
mpm.add_particles_cube(
    cube_origin=[0.7, 0.7, 0.0],
    cube_length=[0.3, 0.3, 0.3],
    material_id=1,
    n_particle_per_cell=4,
    particle_group_id=1
)

# Boundary constraints
mpm.define_boundary_entity()
mpm.add_velocity_constraints(
    [
        {"axis": "x", "bound_loc": "start", "velocity": 0.0},
        {"axis": "x", "bound_loc": "end", "velocity": 0.0},
        {"axis": "y", "bound_loc": "start", "velocity": 0.0},
        {"axis": "y", "bound_loc": "end", "velocity": 0.0},
        {"axis": "z", "bound_loc": "start", "velocity": 0.0},
        {"axis": "z", "bound_loc": "end", "velocity": 0.0}
    ]
)
mpm.add_friction_constrains(
    [
        {"axis": "x", "bound_loc": "start", "sign_n": -1, "friction": 0.4},
        {"axis": "x", "bound_loc": "end", "sign_n": 1, "friction": 0.4},
        {"axis": "y", "bound_loc": "start", "sign_n": -1, "friction": 0.4},
        {"axis": "y", "bound_loc": "end", "sign_n": 1, "friction": 0.4},
        {"axis": "z", "bound_loc": "start", "sign_n": -1, "friction": 0.4},
        {"axis": "z", "bound_loc": "end", "sign_n": 1, "friction": 0.4}
    ]
)

# External loading conditions
mpm.add_external_loadings(
    {"gravity": [0, 0, -9.81]}
)

# Analysis settings
mpm.analysis({
    "mpm_scheme": "usf",
    "locate_particles": False,
    "dt": 1e-05,
    "damping": {
        "type": "Cundall",
        "damping_factor": 0.05
    },
    "resume": {
        "resume": False,
        "step": 0,
        "uuid": "sand3d"
    },
    "velocity_update": False,
    "nsteps": int(1e7),
    "type": "MPMExplicit3D",
    "uuid": "sand3d"
})

# Post-processing
mpm.post_processing({
    "path": "results/",
    "output_steps": 10000,
    "vtk": [
        "displacements"
    ]
})


mpm.write(save_dir=save_dir)

mpm.visualize_mesh(save_path=f'{save_dir}/mesh_config.html', node_indices=True)
mpm.visualize_particles(save_path=f'{save_dir}/particle_config.html')

# Save the current script
# Get the path of the currently running script (main.py)
current_script_path = os.path.abspath(__file__)
utils.save_script(
    current_script_path,
    save_path=f'{save_dir}/input_script.py')
