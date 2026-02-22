
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
