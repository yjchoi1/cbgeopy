import sys
sys.path.append('/work2/08264/baagee/frontera/cbgeopy')
import utils
from mpm import MPMConfig
import demo_utils
from functools import partial
import numpy as np
import vis_utils
import os
import trimesh


save_dir = './'
# Replace the data_dir to the directory containing the obj surface meshes
data_dir = '/work2/08264/baagee/frontera/fundao/fundao_3d-8/'

# Create surface
mesh_base = utils.obj2mesh(f'{data_dir}/mesh_0_v2.obj')  # base
mesh_dike = utils.obj2mesh(f'{data_dir}/mesh_1_v2.obj')  # dike
mesh_sand1 = utils.obj2mesh(f'{data_dir}/mesh_2_v2.obj')  # sand (or slime)
mesh_slime1 = utils.obj2mesh(f'{data_dir}/mesh_2_slime.obj')  # slime-1
mesh_sand_right = utils.obj2mesh(f'{data_dir}/mesh_2_right_sand.obj')  # sand right
mesh_sand2 = utils.obj2mesh(f'{data_dir}/mesh_7.obj')  # sand
mesh_sand3 = utils.obj2mesh(f'{data_dir}/mesh_7_sand.obj')  # sand below slime-2
mesh_slime2 = utils.obj2mesh(f'{data_dir}/mesh_7_slime.obj')  # slime-2
mesh_top = utils.obj2mesh(f'{data_dir}/mesh_8.obj')  # sand
mesh_top_ref = utils.obj2mesh(f'{data_dir}/mesh_top.obj')  # reference top surface for stress initialization


# Make it to trimesh object
floor = demo_utils.gen_surface_mesh(
    [200, 1130], [-850, 450], 10, partial(demo_utils.generate_slope, slope_angle=0, z_min=740, amplitude=0),
    # [180, 1150], [-300, 450], 10, partial(demo_utils.generate_slope, slope_angle=0, z_min=740, amplitude=0),
    plot=False)
mesh_base = trimesh.Trimesh(*mesh_base)
mesh_dike = trimesh.Trimesh(*mesh_dike)
mesh_sand1 = trimesh.Trimesh(*mesh_sand1)
mesh_slime1 = trimesh.Trimesh(*mesh_slime1)
mesh_sand_right = trimesh.Trimesh(*mesh_sand_right)
mesh_sand2 = trimesh.Trimesh(*mesh_sand2)
mesh_sand3 = trimesh.Trimesh(*mesh_sand3)
mesh_slime2 = trimesh.Trimesh(*mesh_slime2)
mesh_top = trimesh.Trimesh(*mesh_top)
mesh_top_ref = trimesh.Trimesh(*mesh_top_ref)


# # See surfaces
# vis_utils.plot_surfaces(floor, mesh_base, mesh_sand_up, mesh_top, save_path='./surface.html')

# Set config
lx, ly, lz = 1130 - 200, 450 - (-850), 950 - 740
mpm = MPMConfig(domain_origin=[200, -850, 740], domain_length=[lx, ly, lz])
# lx, ly, lz = 1150 - 180, 450 - (-300), 950 - 780
# mpm = MPMConfig(domain_origin=[180, -300, 780], domain_length=[lx, ly, lz])

# Mesh
cell_size = 10
mpm.add_mesh(
    n_cells_per_dim=[int(lx/cell_size), int(ly/cell_size), int(lz/cell_size)])

# Add materials
mpm.add_materials(
    [
        {
            "id": 0,
            "name": "liq_sand",
            "density": 2200,
            "youngs_modulus": 4e6,
            "poisson_ratio": 0.30,
            "friction": 0.0,
            "dilation": 0.0,
            "cohesion": 10e3,
            "tension_cutoff": 10,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 45.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb3D"
        },
        {
            "id": 1,
            "density": 2200,
            "name": "slime",
            "youngs_modulus": 2e6,
            "poisson_ratio": 0.30,
            "friction": 0.0,
            "dilation": 0.0,
            "cohesion": 50e3,
            "tension_cutoff": 10,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 45.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb3D"
        },
        {
            "id": 2,
            "name": "bedrock",
            "type": "LinearElastic3D",
            "density": 2600,
            "youngs_modulus": 8e6,
            "poisson_ratio": 0.3
        },
        {
            "id": 3,
            "name": "sand",
            "density": 2200,
            "youngs_modulus": 4e6,
            "poisson_ratio": 0.30,
            "friction": 30.0,
            "dilation": 0.0,
            "cohesion": 1000,
            "tension_cutoff": 100,
            "softening": False,
            "peak_pdstrain": 0.0,
            "residual_friction": 45.0,
            "residual_dilation": 0.0,
            "residual_cohesion": 0.0,
            "residual_pdstrain": 0.0,
            "type": "MohrCoulomb3D"
        },
        {
            "id": 10,
            "name": "geostatic_le",
            "type": "LinearElastic3D",
            "density": 2200,
            "youngs_modulus": 8e6,
            "poisson_ratio": 0.3
        }
    ]
)

n_particle_per_cell = 2
z_find_method = 'linear'
z_fill_method = 'round'
overlap_tolerance = 1
mpm.add_particles_from_topography(
    lower_topography=floor,
    upper_topography=mesh_base,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=0,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_base,
    upper_topography=mesh_dike,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=1,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_dike,
    upper_topography=mesh_sand1,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=2,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_sand1,
    upper_topography=mesh_slime1,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=3,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_slime1,
    upper_topography=mesh_sand_right,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=4,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_sand_right,
    upper_topography=mesh_sand2,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=5,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_sand2,
    upper_topography=mesh_sand3,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=6,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_sand3,
    upper_topography=mesh_slime2,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=7,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
mpm.add_particles_from_topography(
    lower_topography=mesh_slime2,
    upper_topography=mesh_top,
    n_particle_per_cell=n_particle_per_cell,
    material_id=10,
    particle_group_id=8,
    z_find_method=z_find_method,
    base_find_method='simple',
    z_fill_method=z_fill_method,
    overlap_tolerance=overlap_tolerance
)
# mpm.add_particles_from_topography(
#     lower_topography=mesh_top,
#     upper_topography=mesh_top2,
#     n_particle_per_cell=n_particle_per_cell,
#     material_id=10,
#     particle_group_id=8,
#     z_find_method=z_find_method,
#     base_find_method='simple',
#     z_fill_method='round'
# )
mpm.define_particle_entity()

# Stress initialization
mpm.add_initial_stress(
    option='k0',
    top_surface=mesh_top_ref,
    density=2200,
    k0=0.5
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
        {"axis": "x", "bound_loc": "start", "sign_n": -1, "friction": 0.3},
        {"axis": "x", "bound_loc": "end", "sign_n": 1, "friction": 0.3},
        {"axis": "y", "bound_loc": "start", "sign_n": -1, "friction": 0.3},
        {"axis": "y", "bound_loc": "end", "sign_n": 1, "friction": 0.3},
        {"axis": "z", "bound_loc": "start", "sign_n": -1, "friction": 0.3},
        {"axis": "z", "bound_loc": "end", "sign_n": 1, "friction": 0.3}
    ]
)

# External loading conditions
mpm.add_external_loadings(
    {"gravity": [0, 0, -9.81]}
)

# Analysis settings
mpm.analysis({
    "mpm_scheme": "usl",
    "locate_particles": False,
    "dt": 1e-04,
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
    "nsteps": int(1e6),
    "type": "MPMExplicit3D",
    "uuid": "sand3d"
})

# Post-processing
mpm.post_processing({
    "path": "results/",
    "output_steps": 2000,
    "vtk": [
        "displacements",
        "stresses"
    ]
})


mpm.write(save_dir=save_dir)

# mpm.visualize_mesh(save_path=f'{save_dir}/mesh_config.html', node_indices=True)
# mpm.visualize_particles(save_path=f'{save_dir}/particle_config.html')
vis_utils.plot_cross_section(mpm.particle_groups, 'yz', 800, tolerance=3, grid_spacing=cell_size)
# vis_utils.plot_surfaces(floor, mesh_base, mesh_sand_up, mesh_top,
#                         points=mpm.particle_groups,
#                         resolution=1,
#                         save_path=f'{save_dir}/surfaces.html')
vis_utils.save_as_vtk(
    meshes=(floor, mesh_base, mesh_dike, mesh_sand1, mesh_slime1, mesh_sand2, mesh_sand3, mesh_slime2, mesh_top),
    points=mpm.particle_groups
)

# Save the current script
# Get the path of the currently running script (main.py)
current_script_path = os.path.abspath(__file__)
utils.save_script(
    current_script_path,
    save_path=f'{save_dir}/input_script.py')
