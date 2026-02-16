# geometric-pump
https://github.com/tpeng1977/geometric-pump.git

**Suggested GitHub repo name:** `geometric-pump`

Repository for the paper: *Stochastic Thermodynamics and the Falsification of Monotonic Entropy Increase: A Formal Paradox and Empirical Correction* (arXiv).

Nanoscale MD (LAMMPS) simulations and analysis for the "Magic Window" and "Geometric Pump" experiments: isothermal ion systems in asymmetric (funnel) vs symmetric (cylinder) pores to study stochastic entropy and charge separation from thermal fluctuations.

---

## Contents

- **`main.tex`** — LaTeX source (RevTeX 4-2, two-column).
- **`figures/`** — Figures for the manuscript (`time_evolution_statistical.png`, `power_sweep_statistical.png`, `pdf_shift_reconstruction.png`).
- **`script/`** — LAMMPS input generation, run orchestration, and evaluation:
  - `run_magic_window.py` — Multi-seed funnel vs cylinder runs (isothermal NVT, 300 K), optional plotting.
  - `run_power_sweep.py` — Power sweep (retarding potential) runs.
  - `evaluate_magic_window.py` — Compare funnel vs cylinder steady-state net charge and report statistics.
- **`results/`** — Simulation outputs (e.g. `data.funnel_window_s*.dat`, `data.cylinder_control_s*.dat`).
- **`references.bib`** — Bibliography for the paper.

## Simulation setup (LAMMPS)

- **Units:** real (kcal/mol, Å, fs); timestep 0.25 fs.
- **System:** Na⁺ and Cl⁻ in implicit solvent (dielectric 78), LJ + Coulomb; wall/region confinement (cone funnel or cylinder).
- **Ensemble:** Isothermal NVT (Langevin, both reservoirs 300 K).
- **Observables:** Ion counts in left/right track regions → net charge (Na⁺ − Cl⁻) and polarization.

## Quick start

```bash
# Run Magic Window (funnel + cylinder, 4 seeds, 2M steps, 8 workers)
python script/run_magic_window.py

# Evaluate funnel vs cylinder (steady-state bias)
python script/evaluate_magic_window.py

# Compile paper (requires latex + bibtex)
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## License

See repository license file (if present).
