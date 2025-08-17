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

use sang_md::lennard_jones_simulations::{self};

fn main() {
    // main code for running molecular dynamics simulations - version 2
    lennard_jones_simulations::run_md_nve(30, 0.5, 10.0);
}
