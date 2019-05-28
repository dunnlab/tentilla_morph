#!/bin/bash
#SBATCH --job-name=rbmpi_with_constrain

#SBATCH --partition=scavenge
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1

#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00

srun ./revbayes/projects/cmake/rb-mpi revbayes_ct.rev
