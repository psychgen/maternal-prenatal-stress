#!/bin/bash

## Job name (this is what the cluster will show as running - keep short):

#SBATCH --job-name=cortSNPs


#
# Project:

#SBATCH --account=p471


#
# Wall clock limit (change based on resources required):

#SBATCH --time=10:00:00 

#SBATCH --ntasks=1


#
# Output filename customization


#This specification is Jobname_User_JobID

#SBATCH --output=./%x_%u_%j.out


#
# Max memory usage (change based on resources required):

#SBATCH --mem-per-cpu=32G



## Set up job enviroment:
source /cluster/bin/jobsetup
module purge
module load plink/1.90b6.2 
	

# Name used for output file location and filestems (keep short)
	
ANALYSIS_SNAME="cortSNPs" 
	

# Filepath for inputs

INM="/cluster/p/p471/cluster/"
	

# Filepath for outputs
	
OUTM="/cluster/p/p471/cluster/projects/mstress_gxe/"
	
	
#BASE
	
# Amend with filepath/name for your summary statistics file
	
BASE=${OUTM}/canon.txt
	
#TARGET
	
TARGET=${INM}/data/genetic_data/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc

	
#OUTPUT
OUT=${OUTM}

plink --bfile ${TARGET} --double-id  --allow-extra-chr --extract ${BASE} --keep-allele-order --recodeAD include-alt --out ${OUTM}/canon
