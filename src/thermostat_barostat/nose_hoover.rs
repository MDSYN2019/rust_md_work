pub mod nose_hoover {
    use crate::lennard_jones_simulations::{
        compute_pressure_particles, compute_temperature_particles, Particle,
    };

    /// Applies a single-step Nose-Hoover thermostat update to the particle velocities.
    ///
    /// We evolve the thermostat friction variable `xi` with a simple explicit update,
    /// then apply the corresponding exponential damping/scaling to all velocities.
    pub fn apply_thermostat_nose_hoover_particles(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        thermostat_mass: f64,
        dt: f64,
        xi: &mut f64,
    ) {
        if particles.is_empty() || target_temperature <= 0.0 || thermostat_mass <= 0.0 || dt <= 0.0
        {
            return;
        }

        let dof = 3 * particles.len();
        if dof == 0 {
            return;
        }

        let current_temperature = compute_temperature_particles(particles, dof);
        if current_temperature <= 0.0 {
            return;
        }

        // dxi/dt = (T/T0 - 1) * (dof / Q)
        let xi_dot =
            ((current_temperature / target_temperature) - 1.0) * (dof as f64 / thermostat_mass);
        *xi += xi_dot * dt;

        // velocity scaling from friction term v' = v * exp(-xi * dt)
        // clamped for numerical robustness.
        let scale = (-*xi * dt).exp().clamp(0.5, 1.5);
        for p in particles.iter_mut() {
            p.velocity *= scale;
        }
    }

    /// Applies an isotropic Nose-Hoover-like barostat update to particle coordinates
    /// and the simulation box length.
    ///
    /// `eta` is the barostat strain-rate variable (updated in-place).
    pub fn apply_barostat_nose_hoover_particles(
        particles: &mut Vec<Particle>,
        box_length: &mut f64,
        target_pressure: f64,
        barostat_mass: f64,
        dt: f64,
        eta: &mut f64,
    ) {
        if particles.is_empty() || *box_length <= 0.0 || barostat_mass <= 0.0 || dt <= 0.0 {
            return;
        }

        let current_pressure = compute_pressure_particles(particles, *box_length);

        // deta/dt ~ (P - P0) / W
        let eta_dot = (current_pressure - target_pressure) / barostat_mass;
        *eta += eta_dot * dt;

        // isotropic box/coordinate scaling
        let length_scale = (*eta * dt).exp().clamp(0.5, 1.5);
        *box_length *= length_scale;

        for p in particles.iter_mut() {
            p.position *= length_scale;
            // Keep reduced kinetic state consistent with box dilation/compression.
            p.velocity /= length_scale;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::nose_hoover::{
        apply_barostat_nose_hoover_particles, apply_thermostat_nose_hoover_particles,
    };
    use crate::lennard_jones_simulations::{LJParameters, Particle};
    use nalgebra::Vector3;

    fn make_particle(vx: f64, x: f64) -> Particle {
        Particle {
            id: 0,
            position: Vector3::new(x, 0.0, 0.0),
            velocity: Vector3::new(vx, 0.0, 0.0),
            force: Vector3::zeros(),
            lj_parameters: LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 1,
            },
            mass: 1.0,
            energy: 0.0,
            atom_type: 1.0,
            charge: 0.0,
        }
    }

    #[test]
    fn nose_hoover_thermostat_updates_xi_and_velocities() {
        let mut particles = vec![make_particle(5.0, 1.0), make_particle(-5.0, 2.0)];
        let mut xi = 0.0;
        let v_before = particles[0].velocity.norm();

        apply_thermostat_nose_hoover_particles(&mut particles, 1.0, 10.0, 0.01, &mut xi);

        assert!(xi.is_finite());
        assert_ne!(particles[0].velocity.norm(), v_before);
    }

    #[test]
    fn nose_hoover_barostat_scales_box_and_positions() {
        let mut particles = vec![make_particle(1.0, 1.0), make_particle(-1.0, 2.0)];
        let mut box_length = 10.0;
        let mut eta = 0.0;

        apply_barostat_nose_hoover_particles(
            &mut particles,
            &mut box_length,
            0.1,
            5.0,
            0.01,
            &mut eta,
        );

        assert!(box_length > 0.0);
        assert!(eta.is_finite());
        assert!(particles[0].position[0].is_finite());
    }
}
