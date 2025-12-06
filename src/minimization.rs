/*
Energy minimization in molecular dynamics is a preprocessing step
used to find nearby local minimum of potential energy surface

The idea is to remove steric clashes, highly strained geometries,
or unrealistic bond lengths/angbles before running time-evolution simulations
*/

use argmin::prelude::*;
use argmin::solver::quasinewton::BFGS;

pub fn minimize_steepest_descent_particles(
    particles: &mut Vec<Particle>,
    box_length: f64,
    step_size: f64,
    max_iters: usize,
    force_tol: f64,
) {
    use nalgebra::Vector3;

    for iter in 0..max_iters {
        // 1) compute forces from current positions
        compute_forces_particles(particles, box_length);

        // 2) measure max |F|
        let mut max_f = 0.0;
        for p in particles.iter() {
            let f_norm = p.force.norm();
            if f_norm > max_f {
                max_f = f_norm;
            }
        }

        println!("iter {iter:4}   max |F| = {max_f:.6}");

        // 3) check convergence
        if max_f < force_tol {
            println!("Converged: max |F| = {max_f:.6} < {force_tol}");
            break;
        }

        // 4) move along forces: r <- r + α F
        for p in particles.iter_mut() {
            p.position += step_size * p.force;
        }

        // 5) wrap through PBC if needed
        pbc_update(particles, box_length);
    }
}
