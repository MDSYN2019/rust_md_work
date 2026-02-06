pub mod andersen {

    use crate::cell_subdivision;
    use crate::lennard_jones_simulations::{compute_forces_particles, Particle};
    use nalgebra::{zero, Vector3};
    use rand::prelude::*;
    use rand::Rng;
    use rand_distr::{Distribution, Normal};

    pub fn apply_thermostat_andersen_particles(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        dt: f64,
        t_max: f64,
        collision_frequency: f64,
    ) -> () {
        /*
        Initialize system and compute the forces and energy
         */
        let mut t = 0.0;
        let mut switch = 1;

        while target_temperature < t_max {
            // Propagates the half step
            //run_md_andersen_particles(particles, dt, box_length, target_temperature, 1.0, switch);
            //
            //let mut simulation_box = cell_subdivision::SimulationBox {
            //    x_dimension: box_length,
            //    y_dimension: box_length,
            //    z_dimension: box_length,
            //};
            //
            //// Create the subcells - here we have used a subdivision of 10 for the cells
            //let mut subcells = simulation_box.create_subcells(10);
            //// Store the coordinates in cells
            //simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
            //
            //// Compute the forces in the system
            //compute_forces_particles(particles, box_length, &mut subcells);
            //// switches to 2
            //switch = 2;
            //// Propagates the second half time step
            //run_md_andersen_particles(particles, dt, box_length, target_temperature, 1.0, switch);
            //t = t + dt;

            apply_andersen_collisions(particles, target_temperature, collision_frequency, dt);
        }
    }

    pub fn apply_andersen_collisions(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        collision_frequency: f64,
        dt: f64,
    ) -> () {
        /*

        1. Pick a particle

        2. Throw away it's current frequency

        3. Replace it with a new velocity from the maxwell bolztmann distrbiurion
         */
        if dt <= 0.0 || collision_frequency <= 0.0 || target_temperature <= 0.0 {
            return;
        }

        let mut rng = rand::rng();
        let p_coll = 1.0 - (-collision_frequency * dt).exp();

        for p in particles.iter_mut() {
            let r: f64 = rng.random();
            // randomly select a value and change the velocity accordinglty
            if r < p_coll {
                let sigma = (target_temperature / p.mass).sqrt();
                let normal = Normal::new(0.0, sigma).unwrap();

                p.velocity = Vector3::new(
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                );
            }
        }
    }

    pub fn run_md_andersen_particles(
        particles: &mut Vec<Particle>,
        dt: f64,
        _box_length: f64,
        temp: f64,
        nu: f64, // this is the collision frequency
        switch: i64,
    ) -> () {
        // Equations of motion - Andersen thermostat
        let mut a_old: Vec<Vector3<f64>> = Vec::with_capacity(particles.len());
        for a in particles.iter() {
            a_old.push(a.force / a.mass); // compute the acceleration
        }

        if switch == 1 {
            for (p, a_o) in particles.iter_mut().zip(a_old.iter()) {
                // first step velocity verlet
                p.position += dt * p.velocity + (dt * dt) * a_o / 2.0; // update the position current time
                p.velocity += 0.5 * a_o * dt; // update velocity
            }
        } else if switch == 2 {
            /*
            Forces should be recomputed BEFORE this half-kikc
             */
            let mut simulation_box = cell_subdivision::SimulationBox {
                x_dimension: _box_length,
                y_dimension: _box_length,
                z_dimension: _box_length,
            };

            let mut subcells = simulation_box.create_subcells(10);
            // Store the coordinates in cells
            simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
            compute_forces_particles(particles, _box_length, &mut subcells);

            for p in particles.iter_mut() {
                let a_new = p.force / p.mass; // compute the new acceleration
                p.velocity += 0.5 * a_new * dt;
            }
            /*
            Andersen thermostat step - randomize some velocities.

            probability ~ nu * dt ( valid when nu * dt is small; otherwise use 1 - exp(-nu * dt)
             */
            apply_andersen_collisions(particles, temp, nu, dt);
        }
    }
}
