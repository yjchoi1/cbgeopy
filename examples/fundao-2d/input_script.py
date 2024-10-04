import random
import utils
from mpm import MPMConfig
import vis_utils
import os
from tools import random_generator



save_dir = '../layers2d_from_lines/'

# Set config
lx, ly = 1250.0, 150.0
origin_x, origin_y = 0, 0
mpm = MPMConfig(domain_origin=[origin_x, origin_y], domain_length=[lx, ly])

# Mesh
cell_size = 4
mpm.add_mesh(
    n_cells_per_dim=[round(lx/cell_size), round(ly/cell_size)])

# Materials
mpm.add_materials([
    {
        "id": 0,
        "name": "bedrock",
        "type" : "LinearElastic2D",
        "density" : 1800,
        "youngs_modulus" : 50000000,
        "poisson_ratio" : 0.3
    },
    {
        "id": 1,
        "name": "slime",
        "density": 1800,
        "youngs_modulus": 20000000,
        "poisson_ratio": 0.3,
        "friction": 7,
        "dilation": 0.0,
        "cohesion": 1000,
        "tension_cutoff": 100,
        "softening": False,
        "peak_pdstrain": 0.005,
        "residual_friction": 18.0,
        "residual_dilation": 0.0,
        "residual_cohesion": 50.0,
        "residual_pdstrain": 0.05,
        "type": "MohrCoulomb2D"
    },
    {
        "id": 2,
        "name": "sand",
        "density": 1800,
        "youngs_modulus": 20000000,
        "poisson_ratio": 0.3,
        "friction": 10.0,
        "dilation": 0.0,
        "cohesion": 1000,
        "tension_cutoff": 100,
        "softening": False,
        "peak_pdstrain": 0.005,
        "residual_friction": 18.0,
        "residual_dilation": 0.0,
        "residual_cohesion": 50.0,
        "residual_pdstrain": 0.05,
        "type": "MohrCoulomb2D"
    },
    {
        "id": 3,
        "name": "comp_sand",
        "density": 2200,
        "youngs_modulus": 20000000,
        "poisson_ratio": 0.3,
        "friction": 13.0,
        "dilation": 0.0,
        "cohesion": 1000,
        "tension_cutoff": 100,
        "softening": False,
        "peak_pdstrain": 0.005,
        "residual_friction": 18.0,
        "residual_dilation": 0.0,
        "residual_cohesion": 50.0,
        "residual_pdstrain": 0.05,
        "type": "MohrCoulomb2D"
    },
    {
        "id": 100,
        "type": "LinearElastic2D",
        "density": 1800,
        "youngs_modulus": 5e7,
        "poisson_ratio": 0.3
    }
]
)

# Particle
points = [
    [0,	98],
    [37,	97],
    [97,	53],
    [151,	47],
    [423,	40],
    [423,	53],
    [468,	55],
    [491,	40],
    [746,	36],
    [1250,	12],
    [62,	78],
    [410,	78],
    [0,	120],
    [290,	123],
    [354,	102],
    [541,	78],
    [574,	78],
    [655,	56],
    [673,	37],
    [331,	123],
    [380,	104],
    [486,	90],
    [573,	89]
]
mpm.add_particles_from_lines(
    layer_info=[
        {
            "line_points": [points[i] for i in [0, 1, 10, 3, 4, 5, 6, 7, 18, 8, 9]],
            "material_id": 0,
            "particle_group_id": 0
        },
        {
            "line_points": [points[i] for i in [0, 1, 10, 11, 5, 6, 7, 18, 8, 9]],
            "material_id": 1,
            "particle_group_id": 1
        },
        {
            "line_points": [points[i] for i in [12, 13, 14, 15, 16, 17, 18, 8, 9]],
            "material_id": 2,
            "particle_group_id": 2
        },
        {
            "line_points": [points[i] for i in [12, 13, 19, 20, 21, 22, 8, 9]],
            "material_id": 3,
            "particle_group_id": 3
        }
    ],
    n_particle_per_cell=2
)
mpm.remove_overlapping_particles(overlap_tolerance=0.001)
mpm.define_particle_entity()

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
    "nsteps": int(1.8e5),
    "type": "MPMExplicit2D",
    "uuid": "sand2d"
})

# Post-processing
mpm.post_processing({
    "path": "results/",
    "output_steps": 375,
    "vtk": [
        "displacements"
    ]
})

mpm.write(save_dir=save_dir)
vis_utils.plot_scatter(mpm.particle_groups, mpm.domain_ranges, f'{save_dir}/particle_config.png')

# Save the current script
# Get the path of the currently running script (main.py)
current_script_path = os.path.abspath(__file__)
utils.save_script(
    current_script_path,
    save_path=f'{save_dir}/input_script.py')
