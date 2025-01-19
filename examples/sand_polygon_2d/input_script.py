import sys
sys.path.append('/work/08264/baagee/frontera/cbgeopy/')

import random
import utils
from mpm import MPMConfig
import vis_utils
import os
from tools import random_generator
from numpy.random import uniform as randf
from random import randint


# Random parameters
sim_id_range = range(900, 901)

for i in sim_id_range:
    print(f"Generate mpm inputs for simulation {i}...")
    save_dir = f'./sim-{i}'

    # Set config
    lx, ly = 1000.0, 152.0
    origin_x, origin_y = 0, 0
    mpm = MPMConfig(domain_origin=[origin_x, origin_y], domain_length=[lx, ly])

    # Mesh
    cell_size = 4
    mpm.add_mesh(
        n_cells_per_dim=[round(lx/cell_size), round(ly/cell_size)])

    # Materials
    cohesion_options = [10e3, 20e3, 30e3, 40e3, 50e3, 70e3, 100e3, 130e3]
    friction_options = [0, 0, 0, 0, 10, 17.5, 25, 32.5, 40, 43]
    mpm.add_materials([
        {
            "id": 0,
            "type": "LinearElastic2D",
            "density": 1800,
            "youngs_modulus": 50000000.0,
            "poisson_ratio": 0.3
        },
        {
            "id": 1,
            "density": 1800,
            "youngs_modulus": 40e6,
            "poisson_ratio": 0.3,
            "friction": random.choice(friction_options),
            "dilation": 0.0,
            "cohesion": random.choice(cohesion_options),
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        },
        {
            "id": 2,
            "density": 1800,
            "youngs_modulus": 40e6,
            "poisson_ratio": 0.3,
            "friction": random.choice(friction_options),
            "dilation": 0.0,
            "cohesion": random.choice(cohesion_options),
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        },
        {
            "id": 3,
            "density": 1800,
            "youngs_modulus": 40e6,
            "poisson_ratio": 0.3,
            "friction": random.choice(friction_options),
            "dilation": 0.0,
            "cohesion": random.choice(cohesion_options),
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        },
        {
            "id": 4,
            "density": 1800,
            "youngs_modulus": 40e6,
            "poisson_ratio": 0.3,
            "friction": random.choice(friction_options),
            "dilation": 0.0,
            "cohesion": random.choice(cohesion_options),
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 30.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb2D"
        }
    ])

    # Define geometry points
    p1 = [randint(250, 350), 4]
    p2 = [randint(360, 640), 4]
    p3 = [randint(650, 750), 4]
    p5 = [randint(464, 486), randint(12, 70)]
    p7 = [randint(410, 460), randint(74, 85)]
    p8 = [randint(490, 540), randint(74, 95)]

    # Linear interpolation between p7 and p1
    t1 = random.uniform(0.2, 0.8)
    t2 = random.uniform(0.2, 0.8)
    t3 = random.uniform(0.3, 0.7)

    p4 = [p7[i] + t1 * (p1[i] - p7[i]) for i in range(len(p7))]
    p6 = [p8[i] + t2 * (p3[i] - p8[i]) for i in range(len(p8))]
    p78 = [p7[i] + t3 * (p8[i] - p7[i]) for i in range(len(p7))]

    # Particle
    mpm.add_particles_from_lines(
        layer_info=[
            {
                "line_points": [[0, 4], [1000, 4]],
                "material_id": 0,
                "particle_group_id": 0,
                "randomness": 0.0,
            }
        ],
        n_particle_per_cell=2
    )
    mpm.add_particles_from_polygon(
        polygon_info=[
            {
                "polygon_points": [p1, p2, p5, p4],
                "material_id": 1,
                "particle_group_id": 1
            },
            {
                "polygon_points": [p2, p3, p6, p5],
                "material_id": 2,
                "particle_group_id": 2
            },
            {
                "polygon_points": [p4, p5, p78, p7],
                "material_id": 3,
                "particle_group_id": 3
            },
            {
                "polygon_points": [p5, p6, p8, p78],
                "material_id": 4,
                "particle_group_id": 4
            }
        ],
        n_particle_per_cell=2,
        randomness=0.5
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

    vis_utils.plot_scatter(mpm.particle_groups, mpm.domain_ranges, f'{save_dir}/particle_config.png')

    # Save the current script
    # Get the path of the currently running script (main.py)
    current_script_path = os.path.abspath(__file__)
    utils.save_script(
        current_script_path,
        save_path=f'{save_dir}/input_script.py')