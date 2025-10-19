This folder holds **local data** used by the MATLAB pipeline for the paper  
*Task evoked EEG microstate syntax reveals bipolar- and familial- risk-linked network routing beyond static metrics.*

## Paths & configuration

All paths are read from `config/config.m`. Before running the pipeline, edit that file to set:

- `data_root` → `<repo>/data/`
- `data_deriv` → `<repo>/data/derivatives/` (preprocessed & epoched EEG)

- 
## Directory layout expected by the scripts

- **Groups (folder names):** `HC`, `BP_I_Euthymic`, `BP_I_Depressed`, `BP_II_Euthymic`, `BP_II_Depressed`, `Siblings`  
- **Conditions:** `Neutral`, `Negative`, `Positive`  
- **File names:** Any consistent scheme is fine. The examples above are safe:
  - `sub-XXX_task-affect_<Condition>.set` (+ `.fdt` if present)

> If your files live elsewhere or use different names, see the **Manifest option** below.

---

## Manifest option (if you don’t want to rearrange folders)

Create a CSV file (e.g., `data/derivatives_manifest.csv`) like this:


Then adjust the loading block in `code/matlab/Microstate_analysis.m` to iterate over `set_path` instead of scanning directories.  
For portability, prefer **relative** paths or expose the base path in `config/config.m`.

---

## Preprocessing assumptions

The pipeline expects EEG that is:

- **De-identified** and artifact-cleaned (e.g., ICA / artifact rejection completed)
- **Epoched** around the affective task events with consistent condition labels
- Recorded with a consistent channel montage; bad channels interpolated; average reference acceptable
- Sample rate ≥ **250 Hz** recommended

If your preprocessing differs, note it briefly in `docs/METHODS.md`.

---

## What you may commit vs. should not commit

**OK to commit (small, non-identifiable):**
- `participants.tsv` with non-identifiable columns (e.g., `participant_id`, `age`, `sex`, `group`)
- Small **derivative matrices** (e.g., transition matrices, subject-level summary tables) needed to reproduce figures/tables

**Do not commit:**
- Raw EEG (`.eeg`, `.vhdr`, `.vmrk`, `.set`, `.fdt`)  
- Large binary derivatives

> Share large files via a data repository (e.g., OSF or Zenodo) and cite the DOI in the manuscript and `README.md`.

---

## Quick checklist before running

- [ ] `config/config.m` points to your local `data/derivatives` (or a **manifest CSV** is configured)
- [ ] Each subject has three `.set` files (Neutral/Negative/Positive) under the correct **Group/Subject/Condition**
- [ ] Files open cleanly in EEGLAB
- [ ] You have write permissions for `results/` (figures/tables/stats will be created there)

---

## Regenerating results

1. Edit `config/config.m` (paths + toolbox locations).
2. In MATLAB:
   ```matlab
   cd code/matlab
   run_all
