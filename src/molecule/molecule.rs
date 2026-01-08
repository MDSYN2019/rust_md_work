/*

Expanding your rust based molecular dynamics (MD) simulation from point particles
to molecules with bonded interactions and force fields requires adding new types
of interactions

 */

use crate::lennard_jones_simulations::minimum_image_convention;
use crate::lennard_jones_simulations::InitOutput;
use crate::lennard_jones_simulations::LJParameters;
use crate::lennard_jones_simulations::Particle;

use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};

#[derive(Copy, Clone)]
pub struct SimpleBond {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub r0: f64,
}

#[derive(Clone)]
pub struct Atom {
    // Definition of the atom, force on the atom, as well as its current position and velocity
    // as well as the mass and charge
    pub id: usize,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub force: Vector3<f64>,
    pub atom_type: usize,
    pub mass: f64,
    pub charge: f64,
}

#[derive(Clone, Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub k: f64,
    pub r0: f64,
}

pub struct Angle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub k: f64,
}

#[derive(Copy, Clone)]
pub struct NonBondedType {
    pub mass: f64,
    pub charge: f64,
    pub sigma: f64,
    pub epsilon: f64,
}

pub struct ForceField {
    //pub atom_types: HashMap<String, AtomTypeParams>,
    pub bond_types: Vec<Bond>,
    pub angle_types: Vec<Angle>,
}

pub fn bond_distance(a1: &Atom, a2: &Atom) -> f64 {
    // Compute the bond distance between atoms
    (a1.position - a2.position).norm()
}

#[derive(Default, Clone)]
pub struct MoleculeTemplate {
    pub name: String,
    pub atom_types: Vec<String>,      // len
    pub positions: Vec<Vector3<f64>>, // x y z for each atom
    pub bonds: Vec<(usize, usize, f64, f64)>,
    pub exclusion_1_4_scale: Option<f64>, // (i, j, k, k_theta, theta_0)
}

#[derive(Clone, Default, Debug)]
pub struct System {
    pub atoms: Vec<Particle>,
    pub bonds: Vec<Bond>,
}

fn safe_norm(v: &Vector3<f64>) -> f64 {
    let r = v.norm();
    if r < 1e-12 {
        1e-12
    } else {
        r
    }
}

// System is all the atoms (global), bonded terms in global indices, and exclusion sets

pub fn compute_bond_force(atoms: &mut Vec<Particle>, bond: &Bond, box_length: f64) -> f64 {
    /*
    Compute the bond energy,
     */
    let (i, j) = (bond.atom1, bond.atom2); // get atoms#
    let r_vec = atoms[j].position - atoms[i].position; // get vector for position
    let rij_mic = minimum_image_convention(r_vec, box_length);
    let r = rij_mic.norm(); // get distance
    let dr = r - bond.r0; // the difference between the current position and the equilibrium position
    let f_mag = -bond.k * dr; // force magnitude
    let f_vec = (r_vec / r) * f_mag;

    atoms[i].force += f_vec;
    atoms[j].force -= f_vec;

    0.5 * bond.k * dr * dr // return the bond energy
}

pub fn compute_electostatic_bond_short_force(atoms: &mut Vec<Particle>, box_length: f64) -> f64 {
    /*
    Compute the short range real space component of the electrostatic interaction

    https://computecanada.github.io/molmodsim-md-theory-lesson-novice/06-electrostatics/index.html - useful link

    Computing Coulomb potenials is often the most time consuming part of any MD simulation

     */
    let mut total_short_range_potential = 0.0;
    let k_2 = 1.0; // This will be changed to the permittivity of free space
    let e_0 = 1.0;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            // This needs to be properly represent the coloumbing potential - this is a crappy dummy at the moment
            total_short_range_potential += ((atoms[i].charge * atoms[j].charge)
                / (4.0 * 3.14 * e_0))
                / (atoms[0].position - atoms[1].position).norm()
        }
    }
    total_short_range_potential
}

pub fn apply_bonded_forces_and_energy(
    atoms: &mut Vec<Particle>,
    bonds: &[Bond],
    box_length: f64,
) -> f64 {
    /*
    For all the bonds, return the bond energy
     */
    let mut e_bond = 0.0;

    for b in bonds {
        e_bond += compute_bond_force(atoms, b, box_length);
    }
    e_bond
}

pub fn make_h2_system() -> System {
    /*
    Reduced units:

    mass = 1.0 for each H (you can sue 1.0 amu reduced)


    The H-H distance r oscillates around r0 = 0.74

     */
    let r0 = 0.74;
    let k = 100.0;
    let x = 0.5 * r0;

    let mut atoms = vec![
        Particle {
            id: 0,
            position: Vector3::new(-x, 0.0, 0.0),
            velocity: Vector3::new(1.0, 1.0, 1.0),
            force: Vector3::zeros(),
            atom_type: 0.0,
            mass: 1.0,
            charge: 0.0,
            energy: 0.0,
            lj_parameters: (LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 3, // this needs to be corrected
            }),
        },
        Particle {
            id: 1,
            position: Vector3::new(x, 0.0, 0.0),
            velocity: Vector3::new(1.0, -1.0, 1.0),
            force: Vector3::zeros(),
            atom_type: 0.0,
            mass: 1.0,
            charge: 0.0,
            energy: 0.0,
            lj_parameters: (LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 3, // this needs to be corrected
            }),
        },
    ];

    let stretch = 0.05;
    atoms[1].position.x += 0.5 * stretch;
    atoms[0].position.x -= 0.5 * stretch;

    let bonds = vec![Bond {
        atom1: 0,
        atom2: 1,
        k,
        r0,
    }];

    System { atoms, bonds }
}

pub fn create_systems(system: &System, number_of_molecules: i32) -> InitOutput {
    /*
    Create n number of particles
     */
    let mut molecules: Vec<System> = Vec::new();

    for _ in 0..number_of_molecules {
        molecules.push(system.clone());
    }

    // output as the enum we want which will be a valid input to run_md_nve
    InitOutput::Systems(molecules)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_atom() {
        let test_atom = Atom {
            id: 1,
            position: Vector3::new(1.0, 2.0, 3.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            force: Vector3::new(0.0, 0.0, 0.0),
            atom_type: 0,
            mass: 12.01,
            charge: -3.0,
        };
        assert_eq!(test_atom.id, 1);
    }

    #[test]
    fn test_bond_distance_happy() {
        let atom1 = Atom {
            id: 0,
            position: Vector3::new(0.0, 0.0, 0.0),
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            atom_type: 1,
            mass: 1.0,
            charge: 0.0,
        };

        let atom2 = Atom {
            id: 1,
            position: Vector3::new(3.0, 0.0, 0.0),
            ..atom1.clone()
        };

        let dist = bond_distance(&atom1, &atom2);
        assert!((dist - 3.0).abs() < 1e-6) // happy path
    }
}
