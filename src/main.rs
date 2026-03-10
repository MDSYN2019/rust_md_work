//! ----------------------
//! Author: Sang Young Noh
//! ----------------------
//!
//! ------------------------
//! Last Updated: 22/02/2026
//! ------------------------
//!
/*

The HF-self_consistent_field is the standard first-principles
approach for computing approximate quantum mechanical eigenstates
of interacting fermion systems.

Such systems include electrons in atoms, molecules, and condensed matter. Protons
and neutrons in nuclei, and nuclear matter.
*/

use log::error;
#[cfg(feature = "mpi")]
use mpi::traits::*;
use std::env;

use sang_md::lennard_jones_simulations; // this is in lib
use sang_md::molecule::molecule; // this is not in lib - this is the molecule module

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let integrator = env::args()
        .find_map(|arg| arg.strip_prefix("--integrator=").map(str::to_owned))
        .unwrap_or_else(|| "velocity_verlet".to_string());

    let md_mode = match integrator.as_str() {
        "velocity_verlet" => "berendsen",
        "monte_carlo" => "monte_carlo",
        other => {
            log::warn!(
                "Unknown integrator '{other}', falling back to velocity_verlet (Berendsen thermostat)."
            );
            "berendsen"
        }
    };

    // main code for running molecular dynamics simulations - version 2

    #[cfg(feature = "mpi")]
    let universe = mpi::initialize().expect("MPI initialization failed");
    #[cfg(feature = "mpi")]
    let world = universe.world();
    // create a new system
    let mut new_simulation_md =
        match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
            3, 300.0, 30.0, 10.0, 10.0, false,
        ) {
            // How to handle errors - we are returning a result or a string
            Ok(atoms) => atoms,
            Err(e) => {
                error!("Failed to create atoms: {e}");
                return; // Exit early or handle the error as needed
            }
        };
    #[cfg(feature = "mpi")]
    {
        if world.rank() == 0 {
            log::info!("Running MPI-enabled NVE example with integrator={integrator}");
        }
        lennard_jones_simulations::run_md_nve_mpi(
            &mut new_simulation_md,
            30,
            0.0005,
            10.0,
            md_mode,
            &world,
        );
    }

    // When we aren't using MPI distrbuted computing

    #[cfg(not(feature = "mpi"))]
    {
        // running either the default velocity-verlet simulation or monte-carlo simulation
        lennard_jones_simulations::run_md_nve(
            &mut new_simulation_md,
            30,
            0.0005,
            10.0,
            md_mode,
            30.0,
        );
        if md_mode != "monte_carlo" {
            // running an andersen thermostat simulation after velocity-verlet demo
            lennard_jones_simulations::run_md_nve(
                &mut new_simulation_md,
                3000,
                0.0005,
                10.0,
                "andersen",
                30.0,
            );
        }
    }

    // --------------------------------------------------------------------------------------//
    // Create a h2 system
    //  let h2 = molecule::make_h2_system();
    //  let mut systems_vec = molecule::create_systems(&h2, 210);
    //  // assign positions and velocities to the positions
    //
    //  lennard_jones_simulations::set_molecular_positions_and_velocities(&mut systems_vec, 300.0);
    //  // need to modify this - need to implement the create_atoms_with_set_positions_and_velocities to work with molecules here as well
    //  #[cfg(feature = "mpi")]
    //  {
    //      lennard_jones_simulations::run_md_nve_mpi(
    //          &mut systems_vec,
    //          30,
    //          0.0005,
    //          10.0,
    //          "none",
    //          &world,
    //      );
    //  }
    //
    //  #[cfg(not(feature = "mpi"))]
    //  {
    //      lennard_jones_simulations::run_md_nve(&mut systems_vec, 30, 0.0005, 10.0, "none", 3.0);
    //  }
}
