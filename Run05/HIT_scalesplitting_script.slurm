#!/bin/bash

#SBATCH -J HIT_scalesplitting                      # Job name
#SBATCH -o HIT_scalesplitting.o%j                  # Name of stdout output file
#SBATCH -e HIT_scalesplitting.e%j                  # Name of stderr error file
#SBATCH -p skx-dev                            # Queue (partition) name
#SBATCH -N 1                                  # Total # of nodes
#SBATCH --ntasks-per-node 12                              # Total # of cores
#SBATCH -t 02:00:00                           # Run time (hh:mm:ss)
#SBATCH --mail-user=nguyenly@stanford.edu   # Email ID   (UPDATE THIS!)
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A TG-ATM170028                      # Allocation details


# cd to location
cd /work/06005/nguyenly/stampede2/HIT_DNS/Run05
# Launch matlab code...
module load matlab
matlab -nodesktop -nodisplay -nosplash < scalesplitting.m
