/*
Expanding your rust based molecular dynamics (MD) simulation from point particles
to molecules with bonded interactions and force fields requires adding adding new types of interactions
*/
use nalgebra::Vector3;

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
}
