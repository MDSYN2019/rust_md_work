# TODO: Make water-box simulations physically usable for study

## Ensemble correctness (NVT)
- [ ] Replace Berendsen-only production control with an ensemble-correct thermostat for production runs (e.g., properly tuned Nose–Hoover chain or Langevin/BAOAB).
- [ ] Keep Berendsen (or simple velocity-rescaling) only for early equilibration, not production sampling.

## Water model realism
- [ ] Decide target representation explicitly:
  - [ ] Atomistic water (e.g., 3-site/4-site geometry + constraints + electrostatics).
  - [ ] Martini water (correct bead type/parameters, and optional antifreeze/polarizable variants if needed).
- [ ] Add/validate the full nonbonded model required by the chosen representation.

## Electrostatics and long-range interactions
- [ ] Add electrostatic interactions to production-relevant water models.
- [ ] Add long-range electrostatics treatment (PME/Ewald or equivalent) where required.

## Neighbor lists and cutoffs
- [ ] Replace hardcoded cell-pair distance logic with a robust neighbor-list strategy (Verlet list + skin).
- [ ] Validate cutoff, switching/shifting policy, and neighbor-list rebuild frequency.

## Initialization and density control
- [ ] Initialize coordinates using the actual box dimensions (remove hidden fixed-range placement).
- [ ] Build systems at target density from number of molecules/beads and box volume.
- [ ] Remove center-of-mass drift after initialization and periodically during long runs.

## Simulation protocol quality
- [ ] Increase equilibration/production durations to physically meaningful timescales.
- [ ] Add clear run phases: minimization -> equilibration -> production.
- [ ] Record and monitor thermodynamic traces (T, P, energies) for stability.

## Validation and analysis
- [ ] Add radial distribution function (RDF) analysis.
- [ ] Add diffusion coefficient estimation (MSD-based).
- [ ] Validate temperature distribution and fluctuation behavior vs expected ensemble.
- [ ] Compare pressure/density/structure against reference data for the selected model.

## Reproducibility and ergonomics
- [ ] Add seed control and run metadata capture (parameters, code version, date).
- [ ] Add presets for common test systems (small/medium/large water boxes).
- [ ] Document recommended default parameters and known limitations.
