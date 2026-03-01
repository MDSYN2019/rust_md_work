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

#[derive(Clone, Debug)]
pub struct PMEConfig {
    pub alpha: f64,
    pub real_space_cutoff: f64,
    pub grid_dim: usize,
    pub vacuum_permittivity: f64,
}

impl Default for PMEConfig {
    fn default() -> Self {
        Self {
            alpha: 0.35,
            real_space_cutoff: 8.0,
            grid_dim: 16,
            vacuum_permittivity: 1.0,
        }
    }
}

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

#[derive(Clone, Debug)]
pub struct Angle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub k: f64,
    pub theta0: f64,
}

#[derive(Clone, Debug)]
pub struct Dihedral {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub k: f64,
    pub multiplicity: usize,
    pub phase: f64,
}

#[derive(Clone, Debug)]
pub struct Improper {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub k: f64,
    pub psi0: f64,
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
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
    pub impropers: Vec<Improper>,
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

pub fn compute_electostatic_bond_short_force(atoms: &mut Vec<Particle>, _box_length: f64) -> f64 {
    /*
    Compute the short range real space component of the electrostatic interaction

    https://computecanada.github.io/molmodsim-md-theory-lesson-novice/06-electrostatics/index.html - useful link

    Computing Coulomb potenials is often the most time consuming part of any MD simulation

     */
    let mut total_short_range_potential = 0.0;
    let e_0 = 1.0;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let r = (atoms[i].position - atoms[j].position).norm().max(1e-12);
            total_short_range_potential +=
                ((atoms[i].charge * atoms[j].charge) / (4.0 * 3.14 * e_0)) / r
        }
    }
    total_short_range_potential
}

fn pme_flat_index(ix: usize, iy: usize, iz: usize, n: usize) -> usize {
    (ix * n + iy) * n + iz
}

fn erfc_approx(x: f64) -> f64 {
    // Abramowitz and Stegun 7.1.26
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * z);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-(z * z)).exp();
    if x >= 0.0 {
        1.0 - y
    } else {
        1.0 + y
    }
}

fn assign_charge_to_mesh_cic(
    atoms: &[Particle],
    box_length: f64,
    n: usize,
) -> (Vec<f64>, Vec<[usize; 8]>, Vec<[f64; 8]>) {
    let mut rho = vec![0.0; n * n * n];
    let mut particle_nodes = Vec::with_capacity(atoms.len());
    let mut particle_weights = Vec::with_capacity(atoms.len());

    for atom in atoms {
        let fx = atom.position.x.rem_euclid(box_length) / box_length * n as f64;
        let fy = atom.position.y.rem_euclid(box_length) / box_length * n as f64;
        let fz = atom.position.z.rem_euclid(box_length) / box_length * n as f64;

        let ix0 = fx.floor() as usize % n;
        let iy0 = fy.floor() as usize % n;
        let iz0 = fz.floor() as usize % n;
        let ix1 = (ix0 + 1) % n;
        let iy1 = (iy0 + 1) % n;
        let iz1 = (iz0 + 1) % n;

        let tx = fx - ix0 as f64;
        let ty = fy - iy0 as f64;
        let tz = fz - iz0 as f64;

        let nodes = [
            pme_flat_index(ix0, iy0, iz0, n),
            pme_flat_index(ix1, iy0, iz0, n),
            pme_flat_index(ix0, iy1, iz0, n),
            pme_flat_index(ix1, iy1, iz0, n),
            pme_flat_index(ix0, iy0, iz1, n),
            pme_flat_index(ix1, iy0, iz1, n),
            pme_flat_index(ix0, iy1, iz1, n),
            pme_flat_index(ix1, iy1, iz1, n),
        ];

        let weights = [
            (1.0 - tx) * (1.0 - ty) * (1.0 - tz),
            tx * (1.0 - ty) * (1.0 - tz),
            (1.0 - tx) * ty * (1.0 - tz),
            tx * ty * (1.0 - tz),
            (1.0 - tx) * (1.0 - ty) * tz,
            tx * (1.0 - ty) * tz,
            (1.0 - tx) * ty * tz,
            tx * ty * tz,
        ];

        for (node, w) in nodes.iter().zip(weights.iter()) {
            rho[*node] += atom.charge * *w;
        }

        particle_nodes.push(nodes);
        particle_weights.push(weights);
    }

    (rho, particle_nodes, particle_weights)
}

pub fn compute_particle_mesh_ewald(
    atoms: &mut [Particle],
    box_length: f64,
    config: &PMEConfig,
) -> f64 {
    let n = config.grid_dim.max(4);
    let volume = box_length.powi(3);
    let alpha = config.alpha.max(1e-8);
    let k_coulomb = 1.0 / (4.0 * std::f64::consts::PI * config.vacuum_permittivity.max(1e-12));

    let mut total_energy = 0.0;
    let real_cutoff = config.real_space_cutoff.min(0.5 * box_length);

    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let rij = minimum_image_convention(atoms[j].position - atoms[i].position, box_length);
            let r = rij.norm();
            if r <= 1e-12 || r > real_cutoff {
                continue;
            }

            let qij = atoms[i].charge * atoms[j].charge;
            let erfc_term = erfc_approx(alpha * r);
            let prefactor = k_coulomb * qij;
            total_energy += prefactor * erfc_term / r;

            let force_scalar = prefactor
                * (erfc_term / (r * r)
                    + (2.0 * alpha / std::f64::consts::PI.sqrt())
                        * (-(alpha * alpha) * r * r).exp()
                        / r);
            let fij = rij * (force_scalar / r);
            atoms[i].force += fij;
            atoms[j].force -= fij;
        }
    }

    let (rho, particle_nodes, particle_weights) = assign_charge_to_mesh_cic(atoms, box_length, n);

    let mut rho_k_re = vec![0.0; n * n * n];
    let mut rho_k_im = vec![0.0; n * n * n];
    let two_pi = 2.0 * std::f64::consts::PI;

    for kx in 0..n {
        for ky in 0..n {
            for kz in 0..n {
                let mut re = 0.0;
                let mut im = 0.0;
                for ix in 0..n {
                    for iy in 0..n {
                        for iz in 0..n {
                            let idx = pme_flat_index(ix, iy, iz, n);
                            let phase = two_pi * ((kx * ix + ky * iy + kz * iz) as f64 / n as f64);
                            re += rho[idx] * phase.cos();
                            im -= rho[idx] * phase.sin();
                        }
                    }
                }
                let k_idx = pme_flat_index(kx, ky, kz, n);
                rho_k_re[k_idx] = re;
                rho_k_im[k_idx] = im;
            }
        }
    }

    let mut phi_k_re = vec![0.0; n * n * n];
    let mut phi_k_im = vec![0.0; n * n * n];

    for kx in 0..n {
        let kx_i = if kx <= n / 2 {
            kx as i32
        } else {
            kx as i32 - n as i32
        };
        for ky in 0..n {
            let ky_i = if ky <= n / 2 {
                ky as i32
            } else {
                ky as i32 - n as i32
            };
            for kz in 0..n {
                let kz_i = if kz <= n / 2 {
                    kz as i32
                } else {
                    kz as i32 - n as i32
                };
                let k_vec =
                    Vector3::new(kx_i as f64, ky_i as f64, kz_i as f64) * (two_pi / box_length);
                let k_sq = k_vec.norm_squared();
                let idx = pme_flat_index(kx, ky, kz, n);
                if k_sq <= 1e-12 {
                    continue;
                }

                let green = 4.0 * std::f64::consts::PI * (-(k_sq) / (4.0 * alpha * alpha)).exp()
                    / (k_sq * volume);
                phi_k_re[idx] = green * rho_k_re[idx];
                phi_k_im[idx] = green * rho_k_im[idx];

                total_energy += 0.5
                    * k_coulomb
                    * green
                    * (rho_k_re[idx] * rho_k_re[idx] + rho_k_im[idx] * rho_k_im[idx]);
            }
        }
    }

    let mut phi_grid = vec![0.0; n * n * n];
    let norm = 1.0 / (n * n * n) as f64;
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let mut phi = 0.0;
                for kx in 0..n {
                    for ky in 0..n {
                        for kz in 0..n {
                            let idx = pme_flat_index(kx, ky, kz, n);
                            let phase = two_pi * ((kx * ix + ky * iy + kz * iz) as f64 / n as f64);
                            phi += phi_k_re[idx] * phase.cos() - phi_k_im[idx] * phase.sin();
                        }
                    }
                }
                phi_grid[pme_flat_index(ix, iy, iz, n)] = k_coulomb * norm * phi;
            }
        }
    }

    let h = box_length / n as f64;
    let mut ex_grid = vec![0.0; n * n * n];
    let mut ey_grid = vec![0.0; n * n * n];
    let mut ez_grid = vec![0.0; n * n * n];

    for ix in 0..n {
        let ixp = (ix + 1) % n;
        let ixm = (ix + n - 1) % n;
        for iy in 0..n {
            let iyp = (iy + 1) % n;
            let iym = (iy + n - 1) % n;
            for iz in 0..n {
                let izp = (iz + 1) % n;
                let izm = (iz + n - 1) % n;

                let idx = pme_flat_index(ix, iy, iz, n);
                let dphidx = (phi_grid[pme_flat_index(ixp, iy, iz, n)]
                    - phi_grid[pme_flat_index(ixm, iy, iz, n)])
                    / (2.0 * h);
                let dphidy = (phi_grid[pme_flat_index(ix, iyp, iz, n)]
                    - phi_grid[pme_flat_index(ix, iym, iz, n)])
                    / (2.0 * h);
                let dphidz = (phi_grid[pme_flat_index(ix, iy, izp, n)]
                    - phi_grid[pme_flat_index(ix, iy, izm, n)])
                    / (2.0 * h);

                ex_grid[idx] = -dphidx;
                ey_grid[idx] = -dphidy;
                ez_grid[idx] = -dphidz;
            }
        }
    }

    for (atom_idx, atom) in atoms.iter_mut().enumerate() {
        let mut phi_interp = 0.0;
        let mut e_interp = Vector3::zeros();

        for i in 0..8 {
            let node = particle_nodes[atom_idx][i];
            let w = particle_weights[atom_idx][i];
            phi_interp += w * phi_grid[node];
            e_interp.x += w * ex_grid[node];
            e_interp.y += w * ey_grid[node];
            e_interp.z += w * ez_grid[node];
        }

        atom.force += atom.charge * e_interp;
        total_energy += 0.5 * atom.charge * phi_interp;
    }

    let self_energy = atoms.iter().map(|a| a.charge * a.charge).sum::<f64>() * k_coulomb * alpha
        / std::f64::consts::PI.sqrt();

    total_energy - self_energy
}

fn angle_value(atoms: &[Particle], angle: &Angle, box_length: f64) -> f64 {
    let r21 = minimum_image_convention(
        atoms[angle.atom1].position - atoms[angle.atom2].position,
        box_length,
    );
    let r23 = minimum_image_convention(
        atoms[angle.atom3].position - atoms[angle.atom2].position,
        box_length,
    );

    let n1 = r21.norm();
    let n2 = r23.norm();
    if n1 <= 1e-12 || n2 <= 1e-12 {
        return angle.theta0;
    }

    let cos_theta = (r21.dot(&r23) / (n1 * n2)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

fn dihedral_value(atoms: &[Particle], dihedral: &Dihedral, box_length: f64) -> f64 {
    let b1 = minimum_image_convention(
        atoms[dihedral.atom2].position - atoms[dihedral.atom1].position,
        box_length,
    );
    let b2 = minimum_image_convention(
        atoms[dihedral.atom3].position - atoms[dihedral.atom2].position,
        box_length,
    );
    let b3 = minimum_image_convention(
        atoms[dihedral.atom4].position - atoms[dihedral.atom3].position,
        box_length,
    );

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let b2_norm = b2.norm();
    if n1.norm() <= 1e-12 || n2.norm() <= 1e-12 || b2_norm <= 1e-12 {
        return 0.0;
    }

    let b2_hat = b2 / b2_norm;
    let m1 = n1.cross(&b2_hat);

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    y.atan2(x)
}

fn improper_value(atoms: &[Particle], improper: &Improper, box_length: f64) -> f64 {
    let as_dihedral = Dihedral {
        atom1: improper.atom1,
        atom2: improper.atom2,
        atom3: improper.atom3,
        atom4: improper.atom4,
        k: improper.k,
        multiplicity: 1,
        phase: 0.0,
    };
    dihedral_value(atoms, &as_dihedral, box_length)
}

pub fn compute_angle_force(atoms: &mut [Particle], angle: &Angle, box_length: f64) -> f64 {
    let theta = angle_value(atoms, angle, box_length);
    let dtheta = theta - angle.theta0;
    let energy = 0.5 * angle.k * dtheta * dtheta;

    let atom_indices = [angle.atom1, angle.atom2, angle.atom3];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus =
                0.5 * angle.k * (angle_value(atoms, angle, box_length) - angle.theta0).powi(2);
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus =
                0.5 * angle.k * (angle_value(atoms, angle, box_length) - angle.theta0).powi(2);
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn compute_dihedral_force(atoms: &mut [Particle], dihedral: &Dihedral, box_length: f64) -> f64 {
    let phi = dihedral_value(atoms, dihedral, box_length);
    let n = dihedral.multiplicity as f64;
    let energy = dihedral.k * (1.0 + (n * phi - dihedral.phase).cos());

    let atom_indices = [
        dihedral.atom1,
        dihedral.atom2,
        dihedral.atom3,
        dihedral.atom4,
    ];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus = dihedral.k
                * (1.0
                    + ((n * dihedral_value(atoms, dihedral, box_length)) - dihedral.phase).cos());
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus = dihedral.k
                * (1.0
                    + ((n * dihedral_value(atoms, dihedral, box_length)) - dihedral.phase).cos());
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn compute_improper_force(atoms: &mut [Particle], improper: &Improper, box_length: f64) -> f64 {
    let psi = improper_value(atoms, improper, box_length);
    let dpsi = psi - improper.psi0;
    let energy = 0.5 * improper.k * dpsi * dpsi;

    let atom_indices = [
        improper.atom1,
        improper.atom2,
        improper.atom3,
        improper.atom4,
    ];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus = 0.5
                * improper.k
                * (improper_value(atoms, improper, box_length) - improper.psi0).powi(2);
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus = 0.5
                * improper.k
                * (improper_value(atoms, improper, box_length) - improper.psi0).powi(2);
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn apply_all_bonded_forces_and_energy(
    atoms: &mut Vec<Particle>,
    bonds: &[Bond],
    angles: &[Angle],
    dihedrals: &[Dihedral],
    impropers: &[Improper],
    box_length: f64,
) -> f64 {
    let mut energy = 0.0;

    for b in bonds {
        energy += compute_bond_force(atoms, b, box_length);
    }
    for angle in angles {
        energy += compute_angle_force(atoms, angle, box_length);
    }
    for dihedral in dihedrals {
        energy += compute_dihedral_force(atoms, dihedral, box_length);
    }
    for improper in impropers {
        energy += compute_improper_force(atoms, improper, box_length);
    }

    energy
}

pub fn apply_bonded_forces_and_energy(
    atoms: &mut Vec<Particle>,
    bonds: &[Bond],
    box_length: f64,
) -> f64 {
    apply_all_bonded_forces_and_energy(atoms, bonds, &[], &[], &[], box_length)
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

    System {
        atoms,
        bonds,
        angles: vec![],
        dihedrals: vec![],
        impropers: vec![],
    }
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
    #[test]
    fn test_angle_force_energy_zero_at_equilibrium() {
        let mut atoms = vec![
            Particle {
                id: 0,
                position: Vector3::new(1.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 1,
                position: Vector3::new(0.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 2,
                position: Vector3::new(0.0, 1.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
        ];

        let angle = Angle {
            atom1: 0,
            atom2: 1,
            atom3: 2,
            k: 10.0,
            theta0: std::f64::consts::FRAC_PI_2,
        };

        let e = compute_angle_force(&mut atoms, &angle, 10.0);
        assert!(e.abs() < 1e-8);
    }

    #[test]
    fn test_particle_mesh_ewald_two_body_is_finite() {
        let mut atoms = vec![
            Particle {
                id: 0,
                position: Vector3::new(1.0, 1.0, 1.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 1.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 1,
                position: Vector3::new(2.0, 1.0, 1.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: -1.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
        ];

        let config = PMEConfig {
            alpha: 0.4,
            real_space_cutoff: 3.0,
            grid_dim: 4,
            vacuum_permittivity: 1.0,
        };

        let energy = compute_particle_mesh_ewald(&mut atoms, 8.0, &config);
        assert!(energy.is_finite());
        assert!(atoms[0].force.norm().is_finite());
        assert!(atoms[1].force.norm().is_finite());
    }
    #[test]
    fn test_dihedral_energy_phase_shift() {
        let atoms = vec![
            Particle {
                id: 0,
                position: Vector3::new(0.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 1,
                position: Vector3::new(1.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 2,
                position: Vector3::new(1.0, 1.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 3,
                position: Vector3::new(1.0, 1.0, 1.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
        ];

        let dih = Dihedral {
            atom1: 0,
            atom2: 1,
            atom3: 2,
            atom4: 3,
            k: 2.0,
            multiplicity: 1,
            phase: 0.0,
        };

        let phi = dihedral_value(&atoms, &dih, 10.0);
        let e = dih.k * (1.0 + (phi - dih.phase).cos());
        assert!(e.is_finite());
    }
}
