#!/bin/bash
#SBATCH --job-name=revbayes_multitask
#SBATCH --time=120:00:00
#SBATCH --mem=48G
#SBATCH --output=rb_out-%j_.out
#SBATCH --error=rb_out-%j_error.out
#SBATCH --array=1-3
#SBATCH --ntasks-per-node=1
sleep $((SLURM_ARRAY_TASK_ID*60))

filenames = (
	revbayes.rev
	revbayes_ct.rev
	)

file = ${filenames[$SLURM_ARRAY_TASK_ID-1]}
echo $file

./rb $file