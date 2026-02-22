pub mod nose_hoover {
    use crate::lennard_jones_simulations::{
        compute_pressure_particles, compute_temperature_particles,
        Particle,
    };

    pub fn apply_thermostat_nose_hoover_particles(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        thermostat_mass: f64,
        dt: f64,
        xi: &mut f64,
    ) -> () {
        if particles.is_empty() || target_temperature <= 0.0 || thermostat_mass <= 0.0 || dt <= 0.0
        {
            return;
        }

        let dof = 3 * particles.len();
        if dof == 0 {
            return;
        }

        // compute the current temperature

        let current_temperature = compute_temperature_particles(particles, dof);
        if current_temperature <= 0.0 {
            return;
        }
        // dxi/dt = (T/T0 - 1) * (dof / Q)
        let xi_dot =
            ((current_temperature / target_temperature) - 1.0) * (dof as f64 / thermostat_mass);

        *xi += xi_dot * dt;

        // velocity scaling from friction term v' = v * exp(-xi * dt)
        // clamped for numerical robustness

        let scale = (-*xi * dt).exp().clamp(0.5, 1.5);
        for p in particles.iter_mut() {
            p.velocity *= scale;
        }
    }

    /*

    Applies an isotropic Nose-hoover like barostat update to particle coordinates
    and the simulation box length

     */

    pub fn apply_barostat_nose_hoover_particles(
        particles: &mut Vec<Particle>,
        box_length: &mut f64,
        target_pressure: f64,
        barostat_mass: f64,
        dt: f64,
        eta: &mut f64,
    ) -> () {
        if particles.is_empty() || *box_length <= 0.0 || barostat_mass <= 0.0 || dt <= 0.0 {
            return;
        }

        let current_pressure = compute_pressure_particles(particles, *box_length);
        //

        // deta/dt (p-p0)/ W

        let eta_dot = (current_pressure - target_pressure) / barostat_mass;

        *eta += eta_dot * dt;

        // isotropic box/coordinate scaling
        let length_scale = (*eta * dt).exp().clamp(0.5, 1.0);
        *box_length *= length_scale;

        for p in particles.iter_mut() {
            p.position *= length_scale;
            // Keep reduced kinetic state consistent with box dilation/compression
            p.velocity /= length_scale;
        }
    }
}
