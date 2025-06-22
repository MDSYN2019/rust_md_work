//! ----------------------
//! Author: Sang Young Noh
//! ----------------------
//!
//! ------------------------
//! Last Updated: 21/06/2025
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
    //let lj_params_new = lennard_jones_simulations::LJParameters {
    //    epsilon: 1.0,
    //    sigma: 4.0,
    //    number_of_atoms: 2,
    //};
    //
    //let mut new_simulation_md =
    //    match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
    //        3, 10.0, 10.0, 10.0, 10.0,
    //    ) {
    //        // How to handle errors - we are returning a result or a string
    //        Ok(atoms) => atoms,
    //        Err(e) => {
    //            eprintln!("Failed to create atoms: {}", e); // Log the error
    //            return; // Exit early or handle the error as needed
    //        }
    //    };
    //
    //let mut new_simulation_md_clone = new_simulation_md.clone();
    //let mut updated_sim = lennard_jones_simulations::pbc_update(&mut new_simulation_md, 20.0);
    //
    //// Compute the forces
    //lennard_jones_simulations::compute_forces(
    //    &mut new_simulation_md,
    //    lj_params_new.epsilon,
    //    lj_params_new.sigma,
    //);
    //// Running verlet update on the particles
    //lennard_jones_simulations::run_verlet_update(
    //    &mut new_simulation_md,
    //    Vector3::new(0.01, 0.01, 0.01),
    //    0.05,
    //);
    //
    //// compute the temperature of the system
    //let temp = lennard_jones_simulations::compute_temperature(&mut new_simulation_md_clone);
    //println!("The current temperature of the toy system is {}", temp);
    //
    //// apply the thermostat for the system
    //lennard_jones_simulations::apply_thermostat(&mut new_simulation_md_clone, 30.0);
    lennard_jones_simulations::run_md_nve(30, 0.5);
}
