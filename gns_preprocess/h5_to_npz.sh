
#!/bin/bash
module load gcc
module load intel
module load gnuparallel


parallel -j 32 python3 h5_to_npz.py \
--n_dims 2 \
--result_dir "/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers-random/sim-{}/results/sand2d/" \
--mpm_input_path "/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers-random/sim-0/mpm.json" \
--save_path "/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers-random/sim-{}/results/trj.npz" ::: {0..9}
