pub mod berendsen {

    use crate::cell_subdivision;
    use crate::lennard_jones_simulations::compute_pressure_particles; // using the compute_pressure_particles from lib.rs
    use crate::lennard_jones_simulations::Particle; // using the Particle struct from the lennard_jones_simulation mod from lib.rs

    // statistical modules
    use rand::prelude::*;
    use rand::Rng;
    use rand_distr::{Distribution, Normal};

    pub fn apply_barostat_berendsen_particles(
        particles: &mut Vec<Particle>,
        box_length: &mut f64,
        target_pressure: f64,
        tau_p: f64,
        dt: f64,
        compressability: f64,
    ) -> () {
        // if eithe
        if tau_p <= 0.0 || dt <= 0.0 || compressability <= 0.0 || *box_length <= 0.0 {
            return;
        }

        let current_pressure = compute_pressure_particles(particles, *box_length);

        let scale = 1.0 - (dt / tau_p) * compressability * (target_pressure - current_pressure);
        let scale_clamped = scale.clamp(0.5, 1.5); // we bound the scale to be between 0.5 and 1.5
        let length_scale = scale_clamped.cbrt(); // compute the cbrt - this is the final lambda we need to multiple the box vectors and the particle coordinates with

        *box_length *= length_scale;
        // scaling the particle positions
        for p in particles.iter_mut() {
            p.position *= length_scale;
        }
    }
}
