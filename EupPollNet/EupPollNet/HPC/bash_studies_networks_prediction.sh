#! /bin/bash

# these are variables to be used in the job queueing by the HPC:
#$ -q shai.q@bhn27
#$ -cwd
#$ -N embedding_lp_new_evaluation
#$ -l h_vmem=2G
#$ -o ./logs/
#$ -p 1001
#$ -j y

# running the desired R script
/gpfs0/shai/projects/software/R4/R-4.4.1/bin/Rscript HPC_studies_networks_prediction.R $1 
