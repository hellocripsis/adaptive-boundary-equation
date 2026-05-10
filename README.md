# Adaptive Boundary Equation

This repository contains the academic paper, figures, and reproducibility scripts for:

**The Adaptive Boundary Equation: A Minimal Active-Boundary Model for Self-Organizing Information-Entropy Dynamics**

## Contents

paper/    Academic paper PDF and LaTeX source
figures/  Figure files used by the paper
scripts/  Python script used to generate figures

## Reproducing the figures

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/generate_figures.py

The figure-generation script writes figure files into figures/.

## Scope

This repository is for academic paper validation and reproducibility only. It does not contain production Cripsis or Trilithium code, proprietary entropy implementations, customer telemetry, secrets, or operational infrastructure.

## Citation

If you use this work, please cite the paper and this repository.

Archival DOI: https://doi.org/10.5281/zenodo.20090457

## License

Repository code and materials are released under the MIT License unless otherwise noted.
