### Note
# This only works for single material single column collapse

import numpy as np
import pandas as pd
import os
import json


# particles.txt path
path = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/column_geostatic/particles_0.txt'
unit_weight = 1800 * (-9.81)
k0 = 0.5
save_path = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/column_geostatic/'

# Get paritcles
particles = np.genfromtxt(path, skip_header=1)
z_max = particles[:, 2].max()


particle_stress = np.zeros((np.shape(particles)[0], 6))  # second axis is for sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_zx,
particle_stress[:, 0] = k0 * (z_max - particles[:, 2]) * unit_weight  # K0*H*Unit_Weight
particle_stress[:, 1] = k0 * (z_max - particles[:, 2]) * unit_weight  # K0*H*Unit_Weight
particle_stress[:, 2] = (z_max - particles[:, 2]) * unit_weight  # H*Unit_Weight
particle_stress[:, 3] = np.zeros(np.shape(particles)[0])
particle_stress[:, 4] = np.zeros(np.shape(particles)[0])
particle_stress[:, 5] = np.zeros(np.shape(particles)[0])

# Create a header string with the number of particles
header = str(len(particle_stress))

# Save the array to a text file
np.savetxt(f"{save_path}/particles-stresses.txt", particle_stress, delimiter="\t", fmt="%.4f", header=header, comments="")