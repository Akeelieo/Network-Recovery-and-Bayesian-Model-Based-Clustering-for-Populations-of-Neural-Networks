#-------------------------------------------------------------------------------
# make_jobs.R                (repo copy of the Avon working script; verify)
#
# Build a validated SLURM task list (jobs.txt) for nSMC-ABC runs on the
# CHB-MIT data, driven only by patient name(s).
#
# For each seizure it computes a single shared window length L used by BOTH
# the "before" and "during" runs, following the agreed rule:
#
#     L = min( seizure_duration , WINDOW_CAP , available_before )
#
#   where
#     seizure_duration = end - start           (from the summary file)
#     WINDOW_CAP       = 40  seconds            (hard upper cap)
#     available_before = start - prev_end       (clean pre-seizure recording:
#                                                distance from file start, or
#                                                from the previous seizure's end
#                                                in the same file)
#
# Both windows use length L, so before-vs-during is always an equal-length
# comparison for a given seizure. If L < WINDOW_FLOOR the seizure is skipped
# (too little usable recording for meaningful summaries).
#
# Because during-length = L <= seizure_duration, the during-window can never
# run past the end of the file, so no file-length check is needed.
#
# Usage:
#   Rscript make_jobs.R chb02                       # one patient
#   Rscript make_jobs.R chb02 chb05 chb10           # several patients
#   Rscript make_jobs.R chb02 --edf_dir /path --out jobs.txt
#
# Output:
#   jobs.txt        one run per line:  <record> <period> <seizure_index> <L>
#   jobs_report.csv full table of every seizure with its decision (kept/skipped)
#-------------------------------------------------------------------------------
# The summary-parsing function lives in prepare_EEG_functions.R. We only need
# read_seizure_summary(), which has no heavy dependencies, but that file loads
# edf/imputeTS at the top. To avoid requiring those just to build the task list,
# we define a local copy of the parser here (identical logic).
WINDOW_CAP   <- 40   # hard upper cap on window length (seconds)
WINDOW_FLOOR <- 10   # minimum acceptable window length (seconds); below -> skip
#--- argument parsing --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
edf_dir <- "."
out_path <- "jobs.txt"
patients <- character(0)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--edf_dir") { edf_dir <- args[i + 1]; i <- i + 2; next }
  if (a == "--out")     { out_path <- args[i + 1]; i <- i + 2; next }
  patients <- c(patients, a); i <- i + 1
}
if (length(patients) == 0) {
  stop("No patient given. Usage: Rscript make_jobs.R chb02 [chb05 ...] ",
       "[--edf_dir DIR] [--out jobs.txt]")
}
#--- local summary parser (same logic as prepare_EEG_functions.R) ------------
read_seizure_summary <- function(summary_path) {
  if (!file.exists(summary_path)) stop("Summary file not found: ", summary_path)
  lines <- trimws(readLines(summary_path, warn = FALSE))
  current_file <- NA_character_
  starts <- c(); ends <- c(); files <- c()
  for (ln in lines) {
    if (grepl("^File Name:", ln)) {
      current_file <- trimws(sub("^File Name:", "", ln)); next
    }
    if (grepl("^Seizure( [0-9]+)? Start Time:", ln)) {
      v <- as.numeric(gsub("[^0-9.]", "", sub(".*Start Time:", "", ln)))
      starts <- c(starts, v); files <- c(files, current_file); next
    }
    if (grepl("^Seizure( [0-9]+)? End Time:", ln)) {
      v <- as.numeric(gsub("[^0-9.]", "", sub(".*End Time:", "", ln)))
      ends <- c(ends, v); next
    }
  }
  if (length(starts) == 0) stop("No seizures found in ", summary_path)
  if (length(starts) != length(ends))
    stop("Mismatched seizure start/end entries in ", summary_path,
         " (", length(starts), " starts, ", length(ends), " ends)")
  out <- data.frame(file = files,
                    record = sub("\\.edf$", "", files),
                    start = starts, end = ends,
                    duration = ends - starts,
                    stringsAsFactors = FALSE)
  out$seizure_index <- ave(seq_len(nrow(out)), out$file, FUN = seq_along)
  out[, c("file", "record", "seizure_index", "start", "end", "duration")]
}
#--- build the task list -----------------------------------------------------
all_rows <- list()   # per-seizure report rows
job_lines <- character(0)
for (patient in patients) {
  summary_path <- file.path(edf_dir, paste0(patient, "-summary.txt"))
  if (!file.exists(summary_path)) {
    warning("Skipping ", patient, ": summary file not found at ", summary_path,
            call. = FALSE, immediate. = TRUE)
    next
  }
  sz <- read_seizure_summary(summary_path)
  sz <- sz[order(sz$record, sz$start), ]   # ensure chronological within file
  for (k in seq_len(nrow(sz))) {
    rec   <- sz$record[k]
    idx   <- sz$seizure_index[k]
    onset <- sz$start[k]
    dur   <- sz$duration[k]
    # end of the previous seizure in the SAME file (0 if none)
    prev <- sz[sz$record == rec & sz$start < onset, ]
    prev_end <- if (nrow(prev) > 0) max(prev$end) else 0
    available_before <- onset - prev_end
    # the agreed rule
    L <- min(dur, WINDOW_CAP, available_before)
    L <- floor(L)   # whole seconds; keeps grids clean
    # does the edf exist? (only warn; task list can still be built for planning)
    edf_exists <- file.exists(file.path(edf_dir, paste0(rec, ".edf")))
    keep <- is.finite(L) && L >= WINDOW_FLOOR
    reason <- if (!keep) {
      if (!is.finite(L) || L <= 0) "no usable pre-seizure recording"
      else sprintf("L=%.0fs below floor of %ds", L, WINDOW_FLOOR)
    } else if (!edf_exists) {
      "KEPT but .edf missing (download before running)"
    } else "kept"
    all_rows[[length(all_rows) + 1]] <- data.frame(
      patient = patient, record = rec, seizure_index = idx,
      onset = onset, duration = dur,
      available_before = available_before, L = L,
      kept = keep, edf_present = edf_exists, reason = reason,
      stringsAsFactors = FALSE)
    if (keep) {
      # one line per period; both share length L
      job_lines <- c(job_lines,
                     sprintf("%s before %d %d", rec, idx, L),
                     sprintf("%s during %d %d", rec, idx, L))
    }
  }
}
#--- write outputs -----------------------------------------------------------
report <- do.call(rbind, all_rows)
report_path <- paste0(sub("\\.txt$", "", out_path), "_report.csv")
writeLines(job_lines, out_path)
write.csv(report, report_path, row.names = FALSE)
#--- console summary ---------------------------------------------------------
n_seiz  <- nrow(report)
n_kept  <- sum(report$kept)
n_skip  <- n_seiz - n_kept
n_runs  <- length(job_lines)
n_nomiss<- sum(report$kept & !report$edf_present)
cat(sprintf("\nPatients      : %s\n", paste(patients, collapse = ", ")))
cat(sprintf("Seizures found: %d\n", n_seiz))
cat(sprintf("  kept        : %d  -> %d runs (before+during)\n", n_kept, n_runs))
cat(sprintf("  skipped     : %d  (L < %ds floor or no lead-in)\n", n_skip, WINDOW_FLOOR))
if (n_nomiss > 0)
  cat(sprintf("  WARNING     : %d kept seizures have no .edf present yet\n", n_nomiss))
cat(sprintf("\nWrote task list : %s  (%d lines)\n", out_path, n_runs))
cat(sprintf("Wrote report    : %s\n", report_path))
cat(sprintf("\nSubmit with     : sbatch --array=1-%d run_array.sh\n\n", n_runs))
