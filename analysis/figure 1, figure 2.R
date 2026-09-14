# ------------------------------------------------------------------------
# EEG trace figure for the dissertation
#
# Reads chb01_03.edf directly, extracts the four bipolar channels used in
# the analysis, and plots the 40 s before and 40 s during the annotated
# seizure. Writes eeg_trace.pdf for inclusion in the LaTeX document.
# ------------------------------------------------------------------------

## ---- settings ---------------------------------------------------------
edf_path  <- "C:/Users/123/Documents/supp_D_networkABC/chb01_03.edf"
out_pdf   <- "C:/Users/123/Documents/supp_D_networkABC/eeg_trace.pdf"

t_onset   <- 2996          # annotated seizure onset (s), chb01_03
t_offset  <- 3036          # annotated seizure end   (s)
win       <- 40            # seconds shown either side of onset

chan_want <- c("FP1-F7", "FP1-F3", "FP2-F4", "FP2-F8")
scale_fac <- 1             # set to 0.05 to match the model scaling

## ---- package ----------------------------------------------------------
if (!requireNamespace("edfReader", quietly = TRUE)) {
  install.packages("edfReader")
}
library(edfReader)

## ---- read header and match channel names ------------------------------
hdr  <- readEdfHeader(edf_path)
labs <- trimws(hdr$sHeaders$label)

norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
idx  <- match(norm(chan_want), norm(labs))

if (anyNA(idx)) {
  stop("Channel(s) not found: ",
       paste(chan_want[is.na(idx)], collapse = ", "),
       "\nAvailable channels:\n  ", paste(labs, collapse = "\n  "))
}
chan_edf <- labs[idx]          # names exactly as stored in the file

## ---- read the window --------------------------------------------------
t_from <- t_onset - win
t_till <- t_onset + win

sigs <- readEdfSignals(hdr, signals = chan_edf,
                       from = t_from, till = t_till, simplify = FALSE)

Y  <- sapply(chan_edf, function(nm) sigs[[nm]]$signal) * scale_fac
fs <- sigs[[chan_edf[1]]]$sRate

# time axis in seconds relative to seizure onset
tt <- seq(-win, by = 1 / fs, length.out = nrow(Y))

## ---- plot -------------------------------------------------------------
col_left  <- "#1F5FA8"    # matches the 10-20 figure
col_right <- "#B03A2E"
cols      <- c(col_left, col_left, col_right, col_right)

ylim <- range(Y, na.rm = TRUE)

pdf(out_pdf, width = 7.2, height = 5.2, pointsize = 10)
op <- par(mfrow = c(4, 1), mar = c(0.4, 5.0, 0.4, 1.2),
          oma = c(3.8, 0, 1.6, 0), mgp = c(2.8, 0.6, 0), las = 1)

for (k in seq_along(chan_want)) {
  plot(tt, Y[, k], type = "l", col = cols[k], lwd = 0.5,
       xlab = "", ylab = chan_want[k], ylim = ylim,
       xaxt = "n", bty = "n", cex.lab = 1.0, cex.axis = 0.85)
  abline(v = 0, col = "grey30", lty = 2, lwd = 1.2)
  if (k == 1) {
    mtext("before seizure", side = 3, at = -win / 2,
          line = 0.2, cex = 0.75, col = "grey30")
    mtext("during seizure", side = 3, at = win / 2,
          line = 0.2, cex = 0.75, col = "grey30")
  }
  if (k == length(chan_want)) {
    axis(1, cex.axis = 0.85)
    mtext("time relative to seizure onset (s)", side = 1,
          line = 2.4, cex = 0.85)
  }
}
par(op)
dev.off()

cat("written:", out_pdf, "\n")
cat("sampling rate:", fs, "Hz;  samples per channel:", nrow(Y), "\n")