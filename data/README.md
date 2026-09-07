# Data — CHB-MIT Scalp EEG Database

The raw recordings are **not** stored in this repository (tens of gigabytes,
and freely available from PhysioNet). Download them here before running the
pipeline.

Source: CHB-MIT Scalp EEG Database v1.0.0 — https://physionet.org/content/chbmit/1.0.0/
(Shoeb, 2009; distributed via PhysioNet, Goldberger et al., 2000).

## Option A — selective, on-demand (used in this project)

For a given patient, download only its summary and the recordings that contain
a seizure. This is what `hpc/run_patient.sh` does automatically; to do it by
hand for one file:

```
wget -q https://physionet.org/files/chbmit/1.0.0/chb01/chb01-summary.txt
wget -q https://physionet.org/files/chbmit/1.0.0/chb01/chb01_03.edf
```

## Option B — full local mirror

```
# HTTPS
wget -r -N -c -np https://physionet.org/files/chbmit/1.0.0/

# or, anonymously, from the PhysioNet S3 bucket
aws s3 sync --no-sign-request s3://physionet-open/chbmit/1.0.0/ .
```

## Layout expected by the code

The R scripts expect each recording's `.edf` and its patient
`chbNN-summary.txt` to sit in the working directory from which a run is
launched (on the cluster, `~/nJRNMM_edf/`). The summary file supplies the
seizure onset/offset times parsed by `R/prepare_EEG_functions.R`.

Note: `chb24` is present in the recordings but absent from the database's
`SUBJECTINFO` (age/sex) file; `chb17` produced no usable recovered network in
this study.
