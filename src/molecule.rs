/*

Expanding your rust based molecular dynamics (MD) simulation from point particles
to molecules with bonded interactions and force fields requires adding new types
of interactions
 */
use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct SimpleBond {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub r0: f64,
}

#[derive(Clone)]
pub struct Atom {
    pub id: usize,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub force: Vector3<f64>,
    pub atom_type: usize,
    pub mass: f64,
    pub charge: f64,
}

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

pub struct ForceField {
    //pub atom_types: HashMap<String, AtomTypeParams>,
    pub bond_types: Vec<Bond>,
    pub angle_types: Vec<Angle>,
}

pub fn bond_distance(a1: &Atom, a2: &Atom) -> f64 {
    (a1.position - a2.position).norm()
}

fn build_12_exclusions() {}

//pub compute_bond_force(atoms: %mut [Atom], bond: &Bond) {
//
//}

pub fn apply_bonded_forces_and_energy(
    particles: &mut [Particle],
    bonds: &[SimpleBond],
    box_length: f64,
) {
    let (i, j) = (b.i, b.j);
    let r_ij = particles[j].position - particles[i].position; // compute the difference between the positions 
    let r_vec = r_ij;
    let = r_vec.norm();

    if r == 0.0 {continue;}

    let dr = r - b.r0; // difference between the current length and the equilibrium bond length
    e_bond += 0.5 * b.k * dr * dr;
   
    let f_mag = -b.k * dr;
    let f_vec = (r_vec / r) * f_mag;

    particles[i].force  += f_vec;
    particles[j].force += f_vec;
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
