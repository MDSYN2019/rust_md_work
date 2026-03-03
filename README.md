
# rust_md_work  
![CI](https://github.com/MDSYN2019/rust_md_work/actions/workflows/ci.yml/badge.svg)

A high-performance **Molecular Dynamics (MD) engine written in Rust**.  
Implements Lennard-Jones (LJ) interactions, bonded forces, velocity-Verlet time integration, periodic boundary conditions, thermostats, and molecular systems such as H₂.

---

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
