#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=mi_ipw
#
#Project:
#SBATCH --account=p471
#
#Wall clock limit (change based on resources required):
#SBATCH --time=3-00:00:00 
#
#SBATCH --ntasks=1
#
#Output filename customization
#This specification is Jobname_User_JobID
#SBATCH --output=./%x_%u_%j.out
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=50G
#
# Number of processors
#SBATCH --cpus-per-task=6

## Set up job enviroment:

module purge
module load R/4.2.0-foss-2021b


Rscript ./04.3_mi_sib_ipw.R