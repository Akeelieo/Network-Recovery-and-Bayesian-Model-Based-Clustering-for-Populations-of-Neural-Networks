#!/bin/bash
# Repo copy of the working script on Avon (~/nJRNMM_edf/); verify before use.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --output=smc_%x_%j.out
#SBATCH --error=smc_%x_%j.err
#SBATCH --mail-user=u2211111@live.warwick.ac.uk
#SBATCH --mail-type=ALL

. /etc/profile
module purge
module load GCC/13.3.0
module load R/4.4.2
export R_LIBS_USER=~/R/library

cd "$SLURM_SUBMIT_DIR"
Rscript main_SMC_ABC_JRNMM.R $1 $2
