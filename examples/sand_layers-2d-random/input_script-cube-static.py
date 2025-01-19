import random
import utils
from mpm import MPMConfig
import vis_utils
import os
from tools import random_generator
import numpy as np
import copy


# Random parameters
sim_id_range = range(2500, 2800)
n_soil_range = [3, 3]
friction_range = [5, 35]

for i in sim_id_range:
    print(f"Generate mpm inputs for simulation {i}...")
    save_dir = f'./sim-{i}'

    # Set config
    lx, ly = 300.0, 152.0
    origin_x, origin_y = 0, 0
    mpm = MPMConfig(domain_origin=[origin_x, origin_y], domain_length=[lx, ly])

    # Mesh
    cell_size = 4
    mpm.add_mesh(
        n_cells_per_dim=[round(lx/cell_size), round(ly/cell_size)])

    # Materials
    bedrock = [
        {
            "id": 0,
            "type": "LinearElastic2D",
            "density": 1800,
            "youngs_modulus": 50000000.0,
            "poisson_ratio": 0.3
        }
    ]
    soils = random_generator.generate_soils(n_soil_range, friction_range)
    # Add to mpm materials
    materials = bedrock + soils
    mpm.add_materials(materials)

    # Particle
    soil_material_ids = [soil['id'] for soil in soils]
    bedrock_line_points = random_generator.generate_bedrock_line(
        [origin_x, origin_x + lx], [cell_size, 50],
        n_middle_points=random.randint(1, 5)
    )

    x_offset = random.uniform(-50, 50)
    y_offset = np.array(bedrock_line_points).mean(axis=0)[-1]
    mpm.add_particles_from_lines(
        layer_info=[
            {
                "line_points": bedrock_line_points,
                "material_id": 0,
                "particle_group_id": 0,
                "randomness": 0
            },
            {
                "line_points": [
                    [0, 0],
                    [random.uniform(50, 70) + x_offset, 0],
                    [random.uniform(90, 120) + x_offset, random.uniform(30, 60) + y_offset],
                    [random.uniform(180, 210) + x_offset, random.uniform(30, 60) + y_offset],
                    [random.uniform(230, 250) + x_offset, 0],
                    [300, 0]
                ],
                "material_id": random.choice(soil_material_ids),
                "particle_group_id": 1,
                "randomness": 0.3
            },
            {
                "line_points": [
                    [0, 0],
                    [random.uniform(50, 60) + x_offset, 0],
                    [random.uniform(80, 110) + x_offset, random.uniform(50, 90) + y_offset],
                    [random.uniform(190, 220) + x_offset, random.uniform(50, 90) + y_offset],
                    [random.uniform(240, 250) + x_offset, 0],
                    [300, 0]
                ],
                "material_id": random.choice(soil_material_ids),
                "particle_group_id": 2,
                "randomness": 0.3
            }
        ],
        n_particle_per_cell=2
    )
    mpm.remove_overlapping_particles(overlap_tolerance=1)
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
        "nsteps": 180000,
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
