//! ----------------------
//! Author: Sang Young Noh
//! ----------------------
//!
//! ------------------------
//! Last Updated: 02/01/2025
//! ------------------------
//!

/*

The HF-self_consistent_field is the standard first-principles
approach for computing approximate quantum mechanical eigenstates
of interacting fermion systems.

Such systems include electrons in atoms, molecules, and condensed matter. Protons
and neutrons in nuclei, and nuclear matter.
*/

#![allow(unused_variables)] // ensure that unused variables do not cause an error when compiling this program
                            // relax compiler warnings while working through ideas
use std::f64::consts::PI;

use nalgebra::{zero, Vector3};
use sang_md::lennard_jones_simulations;

fn main() {

    // First let's define the force field for the particles in the system 
    let mut lj_params_new = lennard_jones_simulations::LJParameters {
        n: 3,
        i: 0,
        j: 1,
        epsilon: 1.0,
        sigma: 4.0,
        pot: 0.0,
        rij_sq: 0.0,
        sr2: 0.0,
        sr6: 0.0,
        sr12: 0.0,
        epslj: 0.0,
        nsteps: 100,
        na: 2,
    };

    // Running a sample simulation - first generate the velocities and positions of the atoms
    // along with generating 3 molecules within
    let mut new_simulation_md =
        match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
            3, 10.0, 10.0, 10.0,
        ) {
            Ok(atoms) => atoms,
            Err(e) => {
                eprintln!("Failed to create atoms: {}", e); // Log the error
                return; // Exit early or handle the error as needed
            }
        };

    let mut new_simulation_md_clone = new_simulation_md.clone();
    
    // compute the forces - forces are also required to get the acceleration
    lennard_jones_simulations::compute_forces(
        &mut new_simulation_md,
        lj_params_new.epsilon,
        lj_params_new.sigma,
    );

    // running verlet update on the particles
    lennard_jones_simulations::run_verlet_update(
        &mut new_simulation_md_clone,
        Vector3::new(0.01, 0.01, 0.01),
        0.05,
    );

    // compute the temperature of the system
    let temp = lennard_jones_simulations::compute_temperature(&mut new_simulation_md_clone);
    println!("The current temperature of the toy system is {}", temp);

    // apply the thermostat for the system
    lennard_jones_simulations::apply_thermostat(&mut new_simulation_md_clone, 30.0);
}
