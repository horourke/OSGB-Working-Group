#!/bin/sh
#SBATCH --job-name=var0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=josesa@ucr.edu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=2-00:00:00
#SBATCH --output=100_var/experiments1/logs/output%a.out
#SBATCH --array=1-9

pwd; hostname; date

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo This is task $SLURM_ARRAY_TASK_ID
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo 1 Importing R and Rscript:
module load R/4.4.2
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo 3 Running 100_var/103_SimulationScript.R with input $SLURM_ARRAY_TASK_ID 
Rscript 100_var/103_SimulationScript.R $SLURM_ARRAY_TASK_ID 
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

date

