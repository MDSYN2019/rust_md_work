
![image info](FerrumMD.png)


A high-performance **Molecular Dynamics (MD) engine written in Rust**.

Implements Lennard-Jones (LJ) interactions, bonded forces, velocity-Verlet time integration, periodic boundary conditions, thermostats, and molecular systems such as H₂.

---

# FerrumMD
![CI](https://github.com/MDSYN2019/rust_md_work/actions/workflows/ci.yml/badge.svg)

## 🔥 Features

### ✅ Core MD Functionality
- Velocity Verlet integrator  
- Periodic Boundary Conditions (PBC)  
- Minimum Image Convention  
- Site–site Lennard-Jones interactions  
- Bonded interactions via harmonic springs  
- Support for both:
  - **Particle collections** (`InitOutput::Particles`)
  - **Molecular systems** (`InitOutput::Systems`)  

### 🌡 Thermostat Algorithms
- Maxwell–Boltzmann initial velocity sampling  
- Berendsen thermostat (velocity rescaling)  
- Nose-Hoover thermostat (extended-systems coupling)
- Nose-Hoover isotropic barostat (volume/position scaling)
- NVE and pseudo-NVT control  
- Temperature calculation from kinetic energy  

### 🧬 Molecular Support
- Construction of small molecules (e.g., H₂)  
- Support for multiple molecules via system cloning  
- Bonded forces with equilibrium distances and spring constants  

### 🧰 Utilities
- Energy reporting (kinetic, potential, total)  
- Force calculations (bonded + nonbonded LJ)  
- pbc wrapping  
- Configurable time-step, LJ parameters, masses, and box sizes  
- PDB and GRO coordinate readers (`molecule::io::{read_pdb, read_gro}`)  
- Martini `.itp` force-field reader + converter (`molecule::martini::MartiniForceField`)  

---

## 🚀 Getting Started

### Install Rust
You’ll need a stable Rust toolchain:

```bash
rustup update


## ⚡ MPI Parallel NVE Example

An MPI-enabled NVE integration path is available behind the `mpi` feature flag.
It parallelizes Lennard-Jones force/energy accumulation across ranks and uses collective reductions to build global forces.

```bash
cargo run --features mpi
mpirun -n 4 cargo run --features mpi
```

The MPI code path is intended as a parallel-programming example (`run_md_nve_mpi` / `run_md_nve_particles_mpi`).


### Make targets (serial + MPI)

```bash
make run
make run-mpi NP=4
```

## 🧪 Point-particle water-box style trajectory output (GRO + XTC)

This repository now includes a runnable example binary that creates a simple Lennard-Jones point-particle fluid and writes outputs that can be opened in VMD:

```bash
cargo run --bin water_box
```

Generated files:
- `water_box.gro` (final frame structure)
- `water_box.xtc` (trajectory)

### Open in VMD
1. `vmd water_box.gro`
2. In the VMD GUI: **File → Load Data Into Molecule...**
3. Select `water_box.xtc` and load.

> Notes:
> - This is a coarse-grained, point-particle fluid setup (water-like in mass/density intent, not explicit 3-site/4-site water geometry).
> - Coordinates are written in GRO/XTC-compatible units (nm in files).

## 🐍 Python interface (buildable scaffold)

A minimal Python extension interface is available behind the `python` feature.
It exposes:

- `PyMdEngine(sigma, epsilon)` class
- `PyMdEngine.force_at_distance(r)`
- `lj_force_scalar(r, sigma, epsilon)`
- `python_api_version()`

### Build with maturin

```bash
pip install maturin
maturin develop --features python
```

Then in Python:

```python
import sang_md_py

engine = sang_md_py.PyMdEngine(1.0, 1.0)
print(engine.force_at_distance(1.2))
print(sang_md_py.lj_force_scalar(1.2, 1.0, 1.0))
```

This is intended as a starting point you can expand with trajectory stepping, system builders, and observables.

## 🧪 Martini coarse-grained water-box NVT example

A dedicated Martini-style coarse-grained water box example is also available. It runs a single-bead solvent in an NVT-like setup using velocity-Verlet integration with a Berendsen thermostat and writes GRO/XTC outputs:

```bash
cargo run --bin martini_water_box
```

Generated files:
- `martini_water_box.gro`
- `martini_water_box.xtc`

This is intended as a lightweight CG solvent demo that you can visualize in VMD with the same loading flow used for `water_box.xtc`.
