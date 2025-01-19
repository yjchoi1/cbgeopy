import sys
sys.path.append('/home/yj/works/cbgeopy')
import trimesh
import utils
from mpm import MPMConfig
import demo_utils
from functools import partial
import numpy as np
import vis_utils
import os
import trimesh


save_dir = './examples/sand_random_field'

# Random parameters
random_params = {
    "phi_mean": 10.0,
    "phi_std": 2.0,
    "len_scale": 1.0
}

# Set config
lx, ly = 200.0, 60.0
mpm = MPMConfig(domain_origin=[0, 0], domain_length=[lx, ly])

# Mesh
cell_size = 4
mpm.add_mesh(
    n_cells_per_dim=[int(lx/cell_size), int(ly/cell_size)])

# Add materials
bedrock = [
        {
            "id": 0,
            "type": "LinearElastic2D",
            "density": 1800,
            "youngs_modulus": 50000000.0,
            "poisson_ratio": 0.3
        }
    ]

# Particles
# Bedrock
mpm.add_particles_from_lines(
    layer_info=[
        {
            "line_points": [[0, 4], [300, 4]],
            "material_id": 0,
            "particle_group_id": 0
        }
    ],
    n_particle_per_cell=2
)
# Slope
mpm.add_particles_from_random_field(
    polygons_params=[
        {
            "polygon_points": [[0, 4], [60, 4], [60, 20], [0, 20]],
            "random_params": {"mean": 10.0, "std": 2.0, "len_scale": 10.0},
        },
        {
            "polygon_points": [[0, 20], [60, 20], [0, 40]],
            "random_params": {"mean": 20.0, "std": 2.0, "len_scale": 2.0}
        }
    ],
    n_particle_per_cell=2
)
mpm.remove_overlapping_particles(overlap_tolerance=0.001)
mpm.define_particle_entity()

# Material
mpm.add_materials(materials=bedrock)
mpm.add_materials(option="random_field", material_type="MohrCoulomb2D")

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
    "nsteps": 225000,
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

# Save mpm json that will be used after stress initialization
mpm.mpm_json["mesh"]["particles_stresses"] = "particles-stresses.txt"
mpm.write(save_dir=save_dir, file_name='mpm-resume.json')

# mpm json for stress initialization
mpm.mpm_json["mesh"].pop("particles_stresses")
mpm.mpm_json["analysis"]["uuid"] = "sand2d-le"
for particle in mpm.mpm_json["particles"]:
    # Set all materials to bedrock
    particle["generator"]["material_id"] = 0
mpm.mpm_json["analysis"]["resume"]["resume"] = False
mpm.mpm_json["analysis"]["resume"]["uuid"] = "sand2d-le"
mpm.mpm_json["analysis"]["nsteps"] = 70001
# Overwrite
mpm.write(save_dir=save_dir, file_name='mpm-le.json')

# vis_utils.plot_scatter(mpm.particle_groups, mpm.domain_ranges, f'{save_dir}/particle_config.png')
vis_utils.plot_random_field_data(
    mpm.particle_groups, mpm.domain_ranges, f'{save_dir}/random_field_data.png')

# Save the current script
# Get the path of the currently running script (main.py)
current_script_path = os.path.abspath(__file__)
utils.save_script(
    current_script_path,
    save_path=f'{save_dir}/input_script.py')
