#!/bin/bash
#SBATCH --job-name=rb_time_with_constrain

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-23:00:00

./revbayes/projects/cmake/rb timetree_rev/mcmc_TimeTree_siphs_ct.Rev
