# Contributing

Thanks for your interest in improving this project! This repository contains MATLAB code for task-evoked EEG microstate analyses (metrics + directed transition syntax).

## How to propose changes

- **Issues:** Use issues for bugs, feature requests, and questions. Please include:
  - MATLAB version, OS
  - EEGLAB and MicrostateLab versions
  - Exact command(s) run (e.g., `run_all`)
  - Error messages and logs

- **Pull requests (PRs):**
  1. Fork the repo and create a branch from `main` (e.g., `feat/plotting-options`).
  2. Run the pipeline on a small test subset to ensure figures/tables regenerate.
  3. Add/adjust docs (`README`, comments) for any user-facing change.
  4. Open a PR with a short description and screenshots of affected figures.

## Code & style

- MATLAB files should start with a header block describing inputs/outputs and dependencies.
- Replace absolute paths with variables defined in `config/config.m`.
- Use clear function names; avoid modifying globals; prefer returning outputs.
- Keep plotting functions deterministic (set RNG where clustering is involved).

## Data

- Do **not** commit raw EEG. Use `data/` locally.
- Derivative outputs (segmentations, transition matrices) may be shared via release artifacts or data repositories (OSF/Zenodo).

## Reproducibility checklist

- `run_all` completes without errors on your machine.
- New tables/figures land in `results/` and match expected formats.
- Update `README` if the workflow or requirements changed.
- Note versions in `config/versions.md` if toolboxes changed.

## Licensing & citation

- Code is GNU GENERAL PUBLIC LICENSE (see `LICENSE`).
- Cite the accompanying manuscript and this repository (see `CITATION.cff`).
