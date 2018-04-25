#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="p2_test"
#SBATCH --output=p2-srun.out
#SBATCH --mail-user=kimlv@ntnu.no
#SBATCH --mail-type=ALL

mpirun -np 8 ./poisson 128 4

uname -a
