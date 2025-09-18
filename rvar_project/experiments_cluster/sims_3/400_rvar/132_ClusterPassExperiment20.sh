#!/bin/sh
#SBATCH --job-name=rvar2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=josesa@ucr.edu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=2-00:00:00
#SBATCH --output=400_rvar/experiments1/logs/output%a.out
#SBATCH --array=330-649

pwd; hostname; date

echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo This is task $SLURM_ARRAY_TASK_ID,
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
echo 3 Running 400_rvar/103_SimulationScript.R with input $SLURM_ARRAY_TASK_ID 2
Rscript 400_rvar/103_SimulationScript.R $SLURM_ARRAY_TASK_ID 2
echo aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
echo a
echo a

date

