#!/bin/bash
#SBATCH -J mpm_slurm
#SBATCH -o mpm_slurm.o%j
#SBATCH -e mpm_slurm.e%j
#SBATCH -A BCS20003
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p development
#SBATCH --mail-type=all
#SBATCH --mail-user=ychoi402@gatech.edu
#SBATCH -t 2:00:00

#module reset
#module load intel
#module load libfabric


# Inputs
script_dir="/work2/08264/baagee/frontera/cbgeopy/"
result_subdir="/results/sand2d/"
mpm_dir="/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers/"

# Update mpm.json to resume to latest checkpoint if a latest h5 file found.
echo "Update mpm input to resume"
python3 "${script_dir}tools/update_mpm_json.py" \
--json_path="${mpm_dir}/mpm.json" \
--result_dir="${mpm_dir}${result_subdir}" \
--update_option="to_first_checkpoint"

echo "Run MPM from resume point..."
ibrun mpm -i "mpm.json" -f "./"