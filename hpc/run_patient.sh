#!/bin/bash
# Repo copy of the working script on Avon (~/nJRNMM_edf/); verify before use.
# ---------------------------------------------------------------------------
# run_patient.sh  --  submit SMC-ABC before/during jobs for one CHB-MIT patient
#
# Usage:
#   ./run_patient.sh chb11              # all seizure recordings for chb11
#   ./run_patient.sh chb11 03           # only recording chb11_03
#   ./run_patient.sh chb11 03 07 14     # only the listed recordings
#
# What it does, per recording that contains a seizure:
#   1. downloads the patient summary (if not already present)
#   2. finds recordings with >=1 seizure
#   3. downloads each needed .edf (skips if already downloaded)
#   4. clears any stale RefData cache for that recording
#   5. submits a "before" and a "during" job, named chbXX_YY_before / _during
#
# The window cap (min(seizure, 40s)) lives in main_SMC_ABC_JRNMM.R, so this
# script does NOT need to know or set window lengths.
# ---------------------------------------------------------------------------
set -u
BASE="https://physionet.org/files/chbmit/1.0.0"
SBATCH_SCRIPT="run_edf.sbatch"
# --- args ---
patient="${1:-}"
if [ -z "$patient" ]; then
  echo "Usage: $0 <patient> [recording numbers...]"
  echo "  e.g. $0 chb11            (all seizure recordings)"
  echo "       $0 chb11 03 07      (only chb11_03 and chb11_07)"
  exit 1
fi
shift || true
wanted_recordings="$*"     # empty = all seizure recordings
summary="${patient}-summary.txt"
# --- sanity: are we in the right folder? ---
if [ ! -f "$SBATCH_SCRIPT" ]; then
  echo "ERROR: $SBATCH_SCRIPT not found. Run this from ~/nJRNMM_edf."
  exit 1
fi
# --- 1. get the summary ---
if [ ! -f "$summary" ]; then
  echo "Downloading $summary ..."
  wget -q "$BASE/$patient/$summary" || { echo "ERROR: could not download $summary"; exit 1; }
fi
# --- 2. find recordings that contain at least one seizure ---
# lines look like:  File Name: chb11_03.edf   ... Number of Seizures in File: 2
seizure_files=$(awk '
  /^File Name:/ { fn=$3 }
  /^Number of Seizures in File:/ { if ($NF+0 > 0) print fn }
' "$summary" | sed 's/\.edf$//')
if [ -z "$seizure_files" ]; then
  echo "No seizure recordings found in $summary"
  exit 1
fi
echo "Seizure recordings for $patient:"
echo "$seizure_files" | sed 's/^/   /'
echo
# --- helper: is this record wanted? ---
is_wanted() {
  local rec="$1"                    # e.g. chb11_03
  [ -z "$wanted_recordings" ] && return 0
  local num="${rec##*_}"            # 03
  for w in $wanted_recordings; do
    # allow "3" or "03"
    if [ "$num" = "$w" ] || [ "$num" = "0$w" ] || [ "$((10#$num))" = "$((10#$w))" ]; then
      return 0
    fi
  done
  return 1
}
submitted=0
for rec in $seizure_files; do
  is_wanted "$rec" || continue
  edf="${rec}.edf"
  # --- 3. download EDF if missing ---
  if [ ! -s "$edf" ]; then
    echo "Downloading $edf ..."
    wget -q "$BASE/$patient/$edf" || { echo "  WARN: could not download $edf, skipping"; continue; }
  fi
  # --- 4. clear any stale extracted data for this record ---
  rm -rf "RefData/${rec}"
  # --- 5. submit before + during ---
  echo "Submitting $rec (before + during) ..."
  sbatch -J "${rec}_before" "$SBATCH_SCRIPT" "$rec" before
  sbatch -J "${rec}_during" "$SBATCH_SCRIPT" "$rec" during
  submitted=$((submitted + 2))
done
echo
echo "Submitted $submitted jobs for $patient."
echo "Check with:  squeue -u \$USER"
