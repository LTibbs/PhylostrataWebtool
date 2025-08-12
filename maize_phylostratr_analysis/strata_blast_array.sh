#!/bin/bash
#SBATCH -A maizegdb
#SBATCH --job-name="strata_blast"   #name of this job
#SBATCH -p ceres          #name of the partition (queue) you are submitting to
#SBATCH --qos=maizegdb          #name of the partition (queue) you are submitting to
#SBATCH --mem=5GB                  # Real memory (RAM) required (MB), 0 is the whole-node memory
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH -t 02:00:00           #time allocated for this job hours:mins:seconds
#SBATCH -o log/%A_%a_stratablast     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH --mail-user=laura.tibbs-cortes@usda.gov   
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=1-48

date                          #optional, prints out timestamp at the start of the job in stdout file

current_dataset=$1

# load and activate conda environment
module load miniconda
source activate conda_envs/phylostratr/

# submit R script to run phylostratr based on blast results
# the array task ID corresponds to the row number in master_proteomes.csv (after skipping the column names)
Rscript phylostratr_main/strata_blast_array.R $SLURM_ARRAY_TASK_ID $current_dataset

#End of file