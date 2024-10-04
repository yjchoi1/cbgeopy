import random
import utils
from mpm import MPMConfig
import vis_utils
import os
from tools import random_generator


# Random parameters
sim_id_range = range(0, 200)
n_soil_range = [2, 2]
friction_range = [10, 45]

for i in sim_id_range:
    print(f"Generate mpm inputs for simulation {i}...")
    save_dir = f'./sim-{i}'

    # Set config
    lx, ly = 300.0, 150.0
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
    mpm.add_particles_from_lines(
        layer_info=[
            {
                "line_points": [
                    [origin_x, random.uniform(4, 50)] if i % 2 == 0 else
                    [origin_x, 4],
                    [random.uniform(60, 100), 4] if i % 2 == 0 else
                    [random.uniform(200, 240), 4],
                    [origin_x + lx, 4] if i % 2 == 0 else
                    [origin_x + lx, random.uniform(4, 50)],
                ],
                "material_id": 0,
                "particle_group_id": 0
            },
            {
                "line_points": [
                    [origin_x, random.uniform(30, 100)] if i % 2 == 0 else
                    [origin_x, cell_size],

                    [random.uniform(10, 80), random.uniform(70, 100)] if i % 2 == 0 else
                    [random.uniform(180, 220), cell_size],

                    [random.uniform(80, 120), cell_size] if i % 2 == 0 else
                    [random.uniform(220, 290), random.uniform(70, 100)],

                    [origin_x + lx, cell_size] if i % 2 == 0 else
                    [300, random.uniform(30, 100)]
                ],
                "material_id": soil_material_ids[0],
                "particle_group_id": 1
            },
            {
                "line_points": [
                    [origin_x, random.uniform(120, 140)] if i % 2 == 0 else
                    [origin_x, cell_size],

                    [random.uniform(10, 80), random.uniform(80, 120)] if i % 2 == 0 else
                    [random.uniform(170, 220), cell_size],

                    [random.uniform(80, 130), 4] if i % 2 == 0 else
                    [random.uniform(220, 290), random.uniform(80, 120)],

                    [origin_x + lx, cell_size] if i % 2 == 0 else
                    [300, random.uniform(120, 140)]
                ],
                "material_id": soil_material_ids[1],
                "particle_group_id": 2
            },
        ],
        n_particle_per_cell=2,
        randomness=0.8
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
    mpm.add_particle_constraints(
        [
            {
                "pset_id": pid,
                "axis": 'x',
                "velocity": round(random.uniform(0, 0), 2)
            }
            for pid in range(1, len(mpm.particle_groups))
        ]
        +
        [
            {
                "pset_id": pid,
                "axis": 'y',
                "velocity": 0.0
            }
            for pid in range(1, len(mpm.particle_groups))
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
    # # Update input json for resuming without particle velocity constraints
    # current_particle_constraints = mpm.mpm_json['mesh']['boundary_conditions']['particles_velocity_constraints']
    # # Only fix the velocities for the bedrock
    # new_constraints = [constraint for constraint in current_particle_constraints if constraint['pset_id'] == 0]
    # mpm.mpm_json['mesh']['boundary_conditions']['particles_velocity_constraints'] = new_constraints
    mpm.mpm_json['mesh']['boundary_conditions'].pop('particles_velocity_constraints')
    mpm.mpm_json['analysis']['resume']['resume'] = True
    # Overwrite
    mpm.write(save_dir=save_dir, file_name='mpm-resume.json')

    vis_utils.plot_scatter(mpm.particle_groups, mpm.domain_ranges, f'{save_dir}/particle_config.png')

    # Save the current script
    # Get the path of the currently running script (main.py)
    current_script_path = os.path.abspath(__file__)
    utils.save_script(
        current_script_path,
        save_path=f'{save_dir}/input_script.py')
