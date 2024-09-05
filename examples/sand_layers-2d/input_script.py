import trimesh
import utils
from mpm import MPMConfig
import demo_utils
from functools import partial
import numpy as np
import vis_utils
import os
import trimesh


save_dir = './'

# Set config
lx, ly = 200.0, 200.0
mpm = MPMConfig(domain_origin=[0, 0], domain_length=[lx, ly])

# Mesh
cell_size = 5
mpm.add_mesh(
    n_cells_per_dim=[int(lx/cell_size), int(ly/cell_size)])

# Add materials
mpm.add_materials(
    [
        {
            "id": 0,
            "name": "soil-0",
            "density": 1800,
            "youngs_modulus": 20000000.0,
            "poisson_ratio": 0.3,
            "friction": 40.0,
            "dilation": 0.0,
            "cohesion": 100,
            "tension_cutoff": 50,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        },
        {
            "id": 1,
            "name": "soil-2",
            "density": 1800,
            "youngs_modulus": 20000000.0,
            "poisson_ratio": 0.3,
            "friction": 15.0,
            "dilation": 0.0,
            "cohesion": 100,
            "tension_cutoff": 50,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        },
        {
            "id": 10,
            "type": "LinearElastic2D",
            "density": 1800,
            "youngs_modulus": 50000000.0,
            "poisson_ratio": 0.3
        }
    ]
)

# Particles
mpm.add_particles_from_lines(
    layer_info=[
        {
            "line_points": [[0, 50], [80, 10], [200, 30]],
            "material_id": 10,
            "particle_group_id": 0
        }
    ],
    n_particle_per_cell=2
)
mpm.add_particles_cube(
    cube_origin=[40, 40],
    cube_length=[60, 60],
    material_id=0,
    n_particle_per_cell=2,
    particle_group_id=1
)
mpm.add_particles_cube(
    cube_origin=[130, 80],
    cube_length=[60, 60],
    material_id=1,
    n_particle_per_cell=2,
    particle_group_id=2
)
mpm.remove_overlapping_particles(overlap_tolerance=0.001)
mpm.define_particle_entity()
vis_utils.plot_scatter(mpm.particle_groups, mpm.domain_ranges, 'particle_config.png')

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
        {"axis": "x", "bound_loc": "start", "sign_n": -1, "friction": 0.38},
        {"axis": "x", "bound_loc": "end", "sign_n": 1, "friction": 0.38},
        {"axis": "y", "bound_loc": "start", "sign_n": -1, "friction": 0.38},
        {"axis": "y", "bound_loc": "end", "sign_n": 1, "friction": 0.38},
        {"axis": "z", "bound_loc": "start", "sign_n": -1, "friction": 0.38},
        {"axis": "z", "bound_loc": "end", "sign_n": 1, "friction": 0.38}
    ]
)
mpm.add_particle_constraints(
    [
        {
            "pset_id": 0,
            "axis": 'x',
            "velocity": 0.0
        },
        {
            "pset_id": 0,
            "axis": 'y',
            "velocity": 0.0
        },
        {
            "pset_id": 1,
            "axis": 'x',
            "velocity": 10.0
        },
        {
            "pset_id": 1,
            "axis": 'y',
            "velocity": 0.0
        },
        {
            "pset_id": 2,
            "axis": 'x',
            "velocity": -10.0
        },
        {
            "pset_id": 2,
            "axis": 'y',
            "velocity": 0.0
        }
    ]
)

# External loading conditions
mpm.add_external_loadings(
    {"gravity": [0, -9.81]}
)

# Analysis settings
mpm.analysis({
    "mpm_scheme": "usf",
    "locate_particles": False,
    "dt": 1e-04,
    "damping": {
        "type": "Cundall",
        "damping_factor": 0.05
    },
    "resume": {
        "resume": False,
        "step": 0,
        "uuid": "sand2d"
    },
    "velocity_update": False,
    "nsteps": int(1.5e5),
    "type": "MPMExplicit2D",
    "uuid": "sand2d"
})

# Post-processing
mpm.post_processing({
    "path": "results/",
    "output_steps": 375,
    "vtk": [
        "displacements",
        "stresses"
    ]
})

mpm.write(save_dir=save_dir)
# # Update input json for resuming without particle velocity constraints
# current_particle_constraints = mpm.mpm_json['mesh']['boundary_conditions']['particles_velocity_constraints']
# # Only fix the velocities for the bedrock
# new_constraints = [constraint for constraint in current_particle_constraints if constraint['pset_id'] == 0]
# mpm.mpm_json['mesh']['boundary_conditions']['particles_velocity_constraints'] = new_constraints
mpm.mpm_json['mesh']['boundary_conditions'].pop('particles_velocity_constraints')
mpm.mpm_json['analysis']['resume']['resume'] = True
# Overwrite
mpm.write(save_dir=save_dir, file_name='mpm-resume.json')

# Save the current script
# Get the path of the currently running script (main.py)
current_script_path = os.path.abspath(__file__)
utils.save_script(
    current_script_path,
    save_path=f'{save_dir}/input_script.py')
