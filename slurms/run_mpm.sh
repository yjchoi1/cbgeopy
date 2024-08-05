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

# line in mpm.json where resume step is defined
resume_step_line=246
# column in mpm.json where resume step is defined (strictly after `:`)
resume_step_column=14

latest_number=$(ls results/sand3d-1/*.h5 | awk -F'[-_.]' '{print $(NF-1)}' | sort -n | tail -1 | sed 's/^0*//')
sed -i "${resume_step_line}s/^\(.\{${resume_step_column}\}\)[0-9]*\(.*\)$/\1$latest_number\2/" mpm-1.json

ibrun mpm -i "mpm.json" -f "./"