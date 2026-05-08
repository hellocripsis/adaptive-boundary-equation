# Adaptive Boundary Equation

This repository contains the manuscript source, figures, and reproducibility scripts for the paper:

**The Adaptive Boundary Equation: A Minimal Active-Boundary Model for Self-Organizing Information-Entropy Dynamics**

The Adaptive Boundary Equation is a minimal active-boundary model for coupled information-entropy dynamics. The repository provides the LaTeX source, final PDF drafts, generated figures, and Python scripts used to reproduce the numerical checks in the paper.

## Repository structure

```text
paper/    LaTeX source and PDF drafts
figures/  Figure files used by the paper
scripts/  Python scripts used to generate validation figures
outputs/  Generated output figures
```

## Requirements

Python 3 with the packages listed in `requirements.txt`:

```text
numpy
matplotlib
```

## Quickstart

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/generate_figures.py
```

## Reproducing the figures

The main figure-generation script is:

```bash
python scripts/generate_figures.py
```

The critical-slowing validation script is also included separately:

```bash
python scripts/figure5_critical_slowing.py
```

Generated figures are written to `outputs/`.

## Scope and limitations

This repository is for paper validation only. It does not contain production Cripsis or Trilithium code, proprietary entropy implementations, customer telemetry, secrets, or operational infrastructure.

The synthetic demonstrations verify analytic claims of the normalized ABE model. They do not validate any real deployed system.

## Citation

If you use this work, please cite the paper and this repository. A DOI will be added after archival release.

## License

Code and repository materials are released under the MIT License unless otherwise noted.

