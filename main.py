from mpm import MPMConfig


mpm = MPMConfig(domain_origin=[0, 0, 0], domain_length=[1.0, 1.0, 1.0])

# Mesh
mpm.add_mesh(cell_size=[0.1, 0.1, 0.1])

# Add materials
mpm.add_materials(
    [
        {
            "id": 0,
            "density": 1800,
            "youngs_modulus": 2000000.0,
            "poisson_ratio": 0.3,
            "friction": 25.0,
            "dilation": 0.0,
            "cohesion": 100,
            "tension_cutoff": 50,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 42.0,
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
# mpm.add_particles_csv()
mpm.add_particles_cube(
    cube_origin=[0, 0, 0],
    cube_length=[0.4, 0.4, 0.4],
    material_id=0,
    n_particle_per_cell=2,
    particle_group_id=0
)
mpm.add_particles_cube(
    cube_origin=[0.7, 0.7, 0.7],
    cube_length=[0.3, 0.3, 0.3],
    material_id=1,
    n_particle_per_cell=2,
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
    "dt": 1e-06,
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
    "nsteps": 100000,
    "type": "MPMExplicit3D",
    "uuid": "sand3d"
})

# Post-processing
mpm.post_processing({
    "path": "results/",
    "output_steps": 2500,
    "vtk": [
      "displacements"
    ]
})


mpm.write(save_dir='outputs')
a=1