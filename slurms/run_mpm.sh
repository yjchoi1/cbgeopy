#!/bin/bash
#SBATCH -J column_p
#SBATCH -o column_p.o%j
#SBATCH -e column_p.e%j
#SBATCH -A BCS20003
#SBATCH -n 32
#SBATCH -N 4
#SBATCH -p normal
#SBATCH --mail-type=all
#SBATCH --mail-user=ychoi402@gatech.edu
#SBATCH -t 2:00:00

#module reset
#module load intel
#module load libfabric

ibrun mpm -i "mpm.json" -f "./"