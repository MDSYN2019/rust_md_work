//! ----------------------
//! Author: Sang Young Noh
//! ----------------------
//!
//! ------------------------
//! Last Updated: 16/08/2025
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

use sang_md::lennard_jones_simulations;
use sang_md::molecule;

fn main() {
    // main code for running molecular dynamics simulations - version 2
    // create a new system
    let mut new_simulation_md =
        match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
            3, 300.0, 30.0, 10.0, 10.0, false,
        ) {
            // How to handle errors - we are returning a result or a string
            Ok(atoms) => atoms,
            Err(e) => {
                eprintln!("Failed to create atoms: {}", e); // Log the error
                return; // Exit early or handle the error as needed
            }
        };

    lennard_jones_simulations::run_md_nve(&mut new_simulation_md, 30, 0.5, 10.0, "None");

    // Create a h2 system
    let mut h2 = molecule::make_h2_system();
    let mut systems = molecule::create_systems(&h2, 2);

    println!("We have the following atoms {:?}", h2.atoms[0]);
    println!("We have the following atoms {:?}", h2.atoms[1]);

    // need to modify this - need to implement the create_atoms_with_set_positions_and_velocities to work with molecules here as well
    lennard_jones_simulations::run_md_nve(&mut systems, 30, 0.5, 10.0, "None");
}
