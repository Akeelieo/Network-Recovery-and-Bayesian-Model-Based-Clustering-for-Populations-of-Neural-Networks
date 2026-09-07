#!/bin/bash
# Repo copy of the working script on Avon (~/nJRNMM_edf/); verify before use.
#-------------------------------------------------------------------------------
# run_array.sh  --  SLURM array job for nSMC-ABC on CHB-MIT EEG data
#
# Each array task reads ONE line of jobs.txt and runs one inference:
#     <record> <period> <seizure_index> <L>
# e.g.  chb02_16 before 1 40
#
# WORKFLOW
#   1) Generate the task list first:
#          module load R/4.4.2
#          Rscript make_jobs.R chb02
#      -> this writes jobs.txt and tells you how many lines N it has.
#   2) Submit this script as an array of size N:
#          sbatch --array=1-N run_array.sh
#      (replace N with the number make_jobs.R reported, or use the helper
#       submit line at the bottom of this file which counts it for you.)
#
# Each task is INDEPENDENT: if one seizure errors, the others are unaffected.
#-------------------------------------------------------------------------------
#===============================================================================
# SLURM DIRECTIVES  --  values marked ### FILL IN ### must be set for Avon.
# Check them on the cluster with:   sinfo -s        (partitions + time limits)
#===============================================================================
#SBATCH --job-name=nSMCABC
#SBATCH --partition=compute              # e.g. from `sinfo -s`; on Avon often a default partition
#SBATCH --time=24:00:00                  # ### CHECK ### walltime per task (HH:MM:SS); must be <= partition limit
#SBATCH --cpus-per-task=16               # ### CHECK ### cores per task; the R code reads this via SLURM_CPUS_PER_TASK
#SBATCH --mem-per-cpu=2G                 # ### CHECK ### memory per core; raise if runs are killed for OOM
#SBATCH --output=logs/%x_%A_%a.out       # stdout: jobname_arrayjobid_taskid.out
#SBATCH --error=logs/%x_%A_%a.err        # stderr
# NOTE: do NOT hard-code #SBATCH --array here. Pass it at submit time so it
#       matches the current jobs.txt length:   sbatch --array=1-N run_array.sh
#===============================================================================
# ENVIRONMENT
#===============================================================================
set -euo pipefail   # stop on error, undefined var, or failed pipe
module purge
module load GCC/13.3.0
module load R/4.4.2                       # confirmed Avon module string
# Personal R library where SplittingJRNMM and CRAN deps were installed
# on the LOGIN node (compute nodes have no internet).
export R_LIBS_USER=$HOME/R/library
# Make sure the log directory exists (SLURM won't create it for you)
mkdir -p logs
#===============================================================================
# PICK THIS TASK'S LINE FROM jobs.txt
#===============================================================================
JOBS_FILE=jobs.txt
if [[ ! -f "$JOBS_FILE" ]]; then
  echo "ERROR: $JOBS_FILE not found. Run 'Rscript make_jobs.R <patient>' first." >&2
  exit 1
fi
# SLURM_ARRAY_TASK_ID is 1-based; read that line, ignoring blank lines.
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBS_FILE")
if [[ -z "$LINE" ]]; then
  echo "ERROR: no line ${SLURM_ARRAY_TASK_ID} in $JOBS_FILE (is --array range too big?)" >&2
  exit 1
fi
echo "=================================================================="
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Host       : $(hostname)"
echo "Cores      : ${SLURM_CPUS_PER_TASK:-unset}"
echo "Job line   : $LINE"
echo "Started    : $(date)"
echo "=================================================================="
# LINE holds:  record period seizure_index L
# Pass it straight through to the R script as its four arguments.
Rscript main_SMC_ABC_JRNMM.R $LINE
echo "=================================================================="
echo "Finished   : $(date)"
echo "=================================================================="
#-------------------------------------------------------------------------------
# CONVENIENCE: submit with the array range counted automatically
#
#   N=$(grep -cve '^[[:space:]]*$' jobs.txt)   # non-blank lines
#   sbatch --array=1-${N} run_array.sh
#
# Monitor:      squeue -u $USER
# Cancel all:   scancel -u $USER               (or scancel <jobid>)
# One task's log lives in logs/nSMCABC_<jobid>_<taskid>.out
#-------------------------------------------------------------------------------
