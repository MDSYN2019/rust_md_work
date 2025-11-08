/*

=========================================================
 Molecular Dynamics Simulation Framework (Rust)
 Based on:
 "A First Encounter with the Hartree-Fock Self-Consistent Field Method"
 https://apt.scitation.org/doi/abs/10.1119/10.0002644?journalCode=ajp
=========================================================

🔧 Particle Model
-----------------
Each `Particle` struct contains:
- Position: Vector3<f64>
- Velocity: Vector3<f64>
- Force:    Vector3<f64>
- Mass:     f64
- Lennard-Jones Parameters: { sigma, epsilon }

📦 Initialization
-----------------
- Positions initialized randomly inside a cubic simulation box.
- Velocities initialized via the Maxwell-Boltzmann distribution:
    > At thermal equilibrium, particle velocities follow MB statistics.
    > This ensures that kinetic energy corresponds to the target temperature.

🌡️ Temperature Calculation
---------------------------
- Temperature derived from kinetic energy via equipartition theorem:
    T = (2/3) * (KE / N)
  Where:
    KE = total kinetic energy
    N  = number of particles

🌍 Periodic Boundary Conditions (PBC)
-------------------------------------
- Implemented as: xi ← xi mod L
- Mimics an infinite system by wrapping particles across boundaries.
- Prevents artificial wall effects and confinement artifacts.

💥 Lennard-Jones Potential
---------------------------
- Pairwise interaction calculated using:
    V(r) = 4ε [ (σ/r)^12 – (σ/r)^6 ]
- Applied Lorentz–Berthelot mixing rules for ε and σ between different species.

🧊 Velocity Rescaling Thermostat
--------------------------------
- Thermostat applied via:
    λ = sqrt(T_target / T_current)
- Rescales all velocities to control system temperature.
- Simple and stable; not strictly canonical (NVT) but effective for equilibration.

=========================================================

📌 Next Implementation Steps
----------------------------
- [ ] Force calculation from potential (∇V)
- [ ] Center-of-mass velocity removal
- [x] Minimum Image Convention (implemented)
- [x] Energy conservation check (NVE tested)
      → still to analyze for NVT/NPT cases
- [ ] Advanced thermostats:
      → Berendsen (smooth control) - Somewhat done
      → Langevin (stochastic, ensemble-correct) -
- [ ] Radial Distribution Function (RDF)

*/
extern crate assert_type_eq;
mod error;
mod lj_parameters;
mod molecule;

// Use when importing the finished minimization modulexo
//use sang_md::lennard_jones_simulations::{self, compute_total_energy_and_print};

pub mod tensors {
    pub fn outer_product<T>(a: &[T], b: &[T], default_value: T) -> Vec<Vec<T>>
    where
        T: std::ops::Mul<Output = T> + Clone,
        /*
        The function takes two slices &[T] as input. Slices are references
        to arrays or vectors, allowing you to work with portion of a collection,
        to return a vector of vectors of type T
         */
    {
        let mut result = vec![vec![default_value.clone(); b.len()]; a.len()];
        for i in 0..a.len() {
            for j in 0..b.len() {
                result[i][j] = a[i].clone() * b[j].clone();
            }
        }
        result
    }
}

pub mod periodic_boundary_conditions {

    /*
    How do we handle periodic boundaries and minimum image convention in a simulation program?

    Define simulation box and write the necessary methods
    to fix the coordiantes when the molecule has periodic boundary issues

     */
    use nalgebra::Vector3;

    pub struct SimulationBox {
        pub x_dimension: f64,
        pub y_dimension: f64,
        pub z_dimension: f64,
    }

    impl SimulationBox {
        fn cell_subdivison(&self, n_cells: i64) -> () {
            /*

                Cell subdivsision provides a mean for organizing the information about atom positions
            into a form that avoids most of the unnecessary work and reduces the computational effort to a
            O(N_m) level.

            Linked lists are used to assocaite atoms with the cells in which they reside at any given instant.

                 */

            let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);
            let cell_size = box_size / (n_cells as f64);
        }
    }
    pub struct MolecularCoordinates {}
}

#[inline]
pub fn lennard_jones_force_scalar(r: f64, sigma: f64, epsilon: f64) -> f64 {
    // F(r) magnitude along r-hat; positive = repulsive
    // d/dr 4ε[(σ/r)^12 - (σ/r)^6]  =>  24ε [2(σ^12/r^13) - (σ^6/r^7)]
    if r <= 0.0 {
        return 0.0;
    }
    let sr = sigma / r;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    24.0 * epsilon * (2.0 * sr12 - sr6) / r
}

#[inline]
fn safe_norm(x: f64) -> f64 {
    if x < 1e-12 {
        1e-12
    } else {
        x
    }
}

pub mod lennard_jones_simulations {

    use super::*;
    use crate::lj_parameters::lennard_jones_potential;
    use error::compute_average_val;
    use nalgebra::{zero, Vector3};
    use rand::prelude::*;
    use rand::Rng;
    use rand_distr::{Distribution, Normal};

    // importing bonds
    use molecule::apply_bonded_forces_and_energy;
    use molecule::make_h2_system;
    use molecule::Bond;
    use molecule::System;

    #[derive(Clone)]
    pub struct LJParameters {
        // lennard jones parameters and the number of atoms that we have of that parameter
        pub epsilon: f64,
        pub sigma: f64,
        pub number_of_atoms: i32,
    }

    #[derive(Clone)]
    pub struct Particle {
        pub id: usize,
        pub position: Vector3<f64>,
        pub velocity: Vector3<f64>,
        pub force: Vector3<f64>,
        pub lj_parameters: LJParameters,
        pub mass: f64,
        pub energy: f64,
        pub atom_type: f64,
        pub charge: f64,
    }

    #[derive(Clone)]
    pub struct SimulationSummary {
        energy: f64,
    }

    pub enum InitOutput {
        Particles(Vec<Particle>), // define a particles system (single point particle)
        Systems(Vec<System>),     // Define actual molecules
    }

    pub enum InitMode {
        Atoms,
        Molecules,
    }

    impl Particle {
        fn distance(&self, other: &Particle) -> f64 {
            // Compute the distance between two particles
            (self.position - other.position).norm()
        }

        pub fn maxwellboltzmannvelocity(&mut self, temp: f64, mass: f64, _v_max: f64) {
            let mut rng = rand::thread_rng();

            // Standard deviation based on MB distribution
            let sigma_mb = (temp / mass).sqrt();

            // Create a normal distribution with mean = 0, std = sigma
            let normal = Normal::new(0.0, sigma_mb).unwrap();

            // Assign each velocity component independently
            self.velocity[0] = normal.sample(&mut rng);
            self.velocity[1] = normal.sample(&mut rng);
            self.velocity[2] = normal.sample(&mut rng);

            println!("Assigned Maxwell-Boltzmann velocity: {:?}", self.velocity);
        }

        fn update_position_verlet(&mut self, dt: f64) -> () {
            /*
            Verlet scheme to change the position
            Use the verlet scheme to change the velocity
            */
            let a = self.force / self.mass;
            self.position += self.velocity * dt + 0.5 * a * dt * dt;
            //self.position[0] += self.velocity[0] * dt + 0.5 * acceleration[0] * dt * dt;
            //self.position[1] += self.velocity[1] * dt + 0.5 * acceleration[1] * dt * dt;
            //self.position[2] += self.velocity[2] * dt + 0.5 * acceleration[2] * dt * dt;
        }

        fn update_velocity_verlet(&mut self, a_new: Vector3<f64>, dt: f64) {
            /*
            Verlet scheme to update the velocity
             */
            //self.velocity[0] += 0.5 * (old_acceleration[0] + new_acceleration[0]) * dt;
            //self.velocity[1] += 0.5 * (old_acceleration[1] + new_acceleration[1]) * dt;
            //self.velocity[2] += 0.5 * (old_acceleration[2] + new_acceleration[2]) * dt;
            self.velocity += 0.5 * a_new * dt;
        }
    }

    pub fn site_site_energy_calculation(particles: &mut Vec<Particle>, box_length: f64) -> f64 {
        /*
        Computing the total Lennard-Jones energy between all distinct pairs of particles in a molecular system,
        using site-site interactions

        The coordinates r_ia of a site a in molecule i are stored in the elements r(:, i, a)

            For example, if we have two diatomic molecules, then we have r_1a (site a of molecule 1) and
            r_2a (site a of molecule 2). Each molecule is a diatomic molecule (for example, O2).

            We already have a set of particles with the lennard jones parameters defined and stored within. Using
            that data, we need to compute the site_site energy

             */
        let mut total_energy = 0.0;

        for i in 0..particles.len() {
            for j in (i + 1)..particles.len() {
                // double loop over all coordinates in the system

                let sigma_i = particles[i].lj_parameters.sigma; // for particle i, get the sigma
                let epsilon_i = particles[i].lj_parameters.epsilon; // for particle i, get the epsilon
                let sigma_j = particles[j].lj_parameters.sigma; // for particle j, get the sigma
                let epsilon_j = particles[j].lj_parameters.epsilon; // for particle j, get the epsilon

                // Using Lorentz-Bethelot mixing rules
                let computed_sigma = (sigma_i + sigma_j) / 2.0;
                let computed_epsilon = (epsilon_i + epsilon_j).sqrt();
                let r_vec = particles[j].position - particles[i].position;
                let r_vec_mic = minimum_image_convention(r_vec, box_length); // TODO - this needs to be fied
                let r = r_vec_mic.norm();
                let potential = lennard_jones_potential(r, computed_sigma, computed_epsilon);

                // Sum the total energy with the pairwise potential in the system
                total_energy += potential;
            }
        }

        total_energy
    }

    pub fn create_atoms_with_set_positions_and_velocities(
        number_of_atoms: i64,
        temp: f64,
        mass: f64,
        v_max: f64,
        box_dim_max: f64,
        use_atom: bool,
    ) -> Result<InitOutput, String> {
        /*

        Create N atoms with temperature and mass

        Here, we are going with the assumption that we are creating a simulation box
        that is cubic, meaning that we will get a (min + max) * (min + max) *  (min+max) volume
        system for all the molecules

         */
        let mut vector_positions: Vec<Particle> = Vec::new();
        let mut vector_system_positions: Vec<System> = Vec::new();
        let mut rng = rand::rng();
        // Create the number of atoms in the system with the system as necessary
        if !use_atom {
            for index in 0..number_of_atoms {
                // For each particle, we wish to create the initial posiiton,
                // the velocity, and the LJ parameters attached to it
                let mut particle = Particle {
                    // create position for the atom in question
                    position: Vector3::new(
                        // generate x y z position values between -10 and 10
                        rng.random_range(-box_dim_max..box_dim_max),
                        rng.random_range(-box_dim_max..box_dim_max),
                        rng.random_range(-box_dim_max..box_dim_max),
                    ),

                    // create velocity for atom in question
                    velocity: Vector3::new(
                        // generate velocity values between -1 and 1
                        rng.random_range(-1.0..1.0),
                        rng.random_range(-1.0..1.0),
                        rng.random_range(-1.0..1.0),
                    ),
                    // Create the LJ parameters for the atom - default parameters are 1.0, 1.0
                    lj_parameters: (LJParameters {
                        epsilon: 1.0,
                        sigma: 1.0,
                        number_of_atoms: 3,
                    }),
                    force: zero(), // initial force on the atom
                    mass: mass,    // the mass
                    energy: 0.0,
                    atom_type: 0.0,
                    charge: 0.0,
                    id: index as usize,
                };

                // Reset the positions to the maxwell boltzmann distibution of velocities
                particle.maxwellboltzmannvelocity(temp, mass, v_max);
                // push those values into the vector
                vector_positions.push(particle); // push the newly assigned particle into the positions
            }
            Ok(InitOutput::Particles(vector_positions))
        } else {
            // This needs to be fixed
            for _ in 0..number_of_atoms {
                let mut h2_system = make_h2_system(); //
                vector_system_positions.push(h2_system);
            }
            Ok(InitOutput::Systems(vector_system_positions))
        }
    }

    pub fn run_verlet_update_nve(particles: &mut Vec<Particle>, dt: f64, box_length: f64) -> () {
        /*
        Update the position and velocity of the particle using the verlet scheme
         */
        for particle in particles.iter_mut() {
            particle.update_position_verlet(dt);
        }
        pbc_update(particles, box_length);
        compute_forces(particles, box_length);

        for particle in particles.iter_mut() {
            println!(
                "The original position and velocity is {:?} and {:?} ",
                particle.position, particle.velocity
            );
            let a_new = particle.force / particle.mass;
            particle.update_velocity_verlet(a_new, dt);

            println!(
                "After a iteration step, the position and velocity is {:?} and {:?} ",
                particle.position, particle.velocity
            );
        }
    }

    pub fn apply_bond_force(particles: &mut [Particle], b: &Bond, box_length: f64) -> f64 {
        let rij = particles[b.atom1].position - particles[b.atom2].position;
        let rij_mic = minimum_image_convention(rij, box_length);
        let r = safe_norm(rij_mic.norm());
        let dr = r - b.r0;
        let f_mag = -b.k * dr; // along r̂, attractive if r>r0
        let f_vec = (rij_mic / r) * f_mag; // vector force on i

        particles[b.atom1].force += f_vec;
        particles[b.atom2].force -= f_vec;

        0.5 * b.k * dr * dr
    }

    pub fn compute_forces(particles: &mut Vec<Particle>, bonds: &[Bond], box_length: f64) {
        /*
        Computing forces between the molecules
         */
        for p in particles.iter_mut() {
            p.force = Vector3::zeros();
        }

        let n = particles.len(); // number of particles in the system
                                 // initalize zero forces for each particle
        for i in 0..n {
            for j in (i + 1)..n {
                let r_vec = particles[j].position - particles[i].position;
                let r_mic = minimum_image_convention(r_vec, box_length);
                let r = r_mic.norm();
                if r == 0.0 {
                    continue;
                }

                // mix params (Lorentz-Berthelot)
                let si = particles[i].lj_parameters.sigma;
                let ei = particles[i].lj_parameters.epsilon;
                let sj = particles[j].lj_parameters.sigma;
                let ej = particles[j].lj_parameters.epsilon;
                let sigma = 0.5 * (si + sj);
                let epsilon = (ei * ej).sqrt();
                let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
                let f_vec = (r_mic / r) * f_mag; // along r-hat

                // action = -reaction
                particles[i].force -= f_vec;
                particles[j].force += f_vec;
                println!(
                    "The forces are {:?} {:?}",
                    particles[i].force, particles[j].force
                );
            }
        }
        // apply bonded terms
        let _ = apply_bonded_forces_and_energy();
    }

    pub fn compute_temperature(particles: &mut Vec<Particle>, dof: usize) -> f64 {
        /*
            Compute the current temperature of the system

            To implement a thermostat for a molecular dynamics (MD) simulation in Rust, we need
            to control the temperature of the system by adjusting the velocities of the particles

            This is typically done using methods like Berendson thermostat or velocity rescaling

            ---

            T = 2/3 * (KE/N) - this is the fundamental relation in classical statistical mechanics that connects
            the temperature of a system with its kinetic energy per particle

        ---

        dof is used to account of degrees of freedom, to account for center of mass drift of the system
         */
        if dof == 0 {
            return 0.0;
        }
        let mut total_kinetic_energy = 0.0;
        let num_particles = particles.len() as f64;
        for particle in particles {
            let velocity_sq = particle.velocity.norm_squared();
            total_kinetic_energy += 0.5 * particle.mass * velocity_sq;
        }
        // T = (2/3) * (KE / N)
        (2.0) * (total_kinetic_energy) / (dof as f64)
    }

    pub fn apply_thermostat(particles: &mut Vec<Particle>, target_temperature: f64) -> () {
        /*
            Here we are rescaling the velocity in the current step towards the target temperature

        The full term for computing the thermostat is as follows:

        T=3NkB​2​i∑​21​mvi2​


         */
        // computing the degrees of freedom

        let dof = 3 * particles.len().saturating_sub(3);

        let current_temperature = compute_temperature(particles, dof);
        if current_temperature == 0.0 {
            return; // Avoid division by zero
        }

        // Compute scaling factor
        let lambda = (target_temperature / current_temperature).sqrt();

        // Rescale velocities
        for particle in particles {
            particle.velocity *= lambda;
        }

        println!(
            "Applied thermostat: Current Temp = {:.2}, Target Temp = {:.2}, Scaling Factor = {:.2}",
            current_temperature, target_temperature, lambda
        );
    }

    pub fn apply_thermostat_berendsen(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        tau: f64,
        dt: f64,
    ) -> () {
        /*

        If the system's instantanepus temperature T differs from the target temperature T_0, the Berendsen
        thermostat weakly couples the system to a head bath that gently nudges T towards T_0 over a characteristic
        relaxation time

         */
        let dof = 3 * particles.len().saturating_sub(3);
        let current_temperature = compute_temperature(particles, dof);

        //if tau <= 0.0 || dt <= 0.0 || current_temperature <= 0.0 || target_temperature <= 0.0 {
        //    return 1.0;
        //}

        // Discrete Berendsen: T' = T * (1 + (dt/tau)(T_0/T - 1))
        // Velocitiea s scale as sqrt (T'/T)
        let x = (dt / tau) * (target_temperature / current_temperature - 1.0);

        // clamp to avoid negative
        let x_clamped = x.clamp(-0.9, 10.0);
        let lambda = (1.0 + x_clamped).max(1e-12).sqrt();

        for particle in particles {
            particle.velocity *= lambda;
        }
    }

    pub fn pbc_update(particles: &mut Vec<Particle>, box_length: f64) {
        for particle in particles.iter_mut() {
            for i in 0..3 {
                particle.position[i] = particle.position[i].rem_euclid(box_length);
            }
        }
    }

    pub fn compute_total_energy_and_print(particles: &mut Vec<Particle>, box_length: f64) -> f64 {
        /*
        compute the total kinetic + potential energy of the system
         */
        let mut kinetic_energy = 0.0;

        for p in particles.iter_mut() {
            let v2 = p.velocity.norm_squared();
            kinetic_energy += 0.5 * p.mass * v2;
        }

        let potential_energy = site_site_energy_calculation(particles, box_length);

        println!(
            "The potential energy is {}",
            kinetic_energy + potential_energy
        );

        let mut bond_energ = 0.0;

        kinetic_energy + potential_energy
    }

    pub fn minimum_image_convention(rij: Vector3<f64>, box_length: f64) -> Vector3<f64> {
        Vector3::new(
            rij[0] - box_length * (rij[0] / box_length).round(),
            rij[1] - box_length * (rij[1] / box_length).round(),
            rij[2] - box_length * (rij[2] / box_length).round(),
        )
    }

    pub fn run_md_nve(
        system: &mut Vec<Particle>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
    ) {
        /*
        We are now equipt to implement a NVE molecular dynamics simulations.
        define time step and number of steps
         */
        let mut final_summary = SimulationSummary { energy: 0.0 };

        let lj_params_new = LJParameters {
            epsilon: 1.0,
            sigma: 4.0,
            number_of_atoms: 2,
        };

        let mut values: Vec<f32> = Vec::new();

        // Compute the initial total energy of the system
        let initial_energy = compute_total_energy_and_print(system, box_length);
        // Loop over the total system for number_of_steps
        for i in 0..number_of_steps {
            // update periodic boundary conditions
            pbc_update(system, box_length);
            // update velocities using the verlet format
            run_verlet_update_nve(system, 0.05, box_length);
            let dof = 3 * system.len().saturating_sub(3);
            let system_temperature = compute_temperature(system, dof);
            println!("The temperature of the system is {}", system_temperature);

            if thermostat == "berendsen" {
                apply_thermostat_berendsen(system, 300.0, 0.1, 0.05);
            } else {
                apply_thermostat(system, 300.0);
            }
            let total_energy = compute_total_energy_and_print(system, box_length);
            // update the summary

            final_summary.energy = total_energy;
            values.push(total_energy as f32);
        }

        compute_average_val(&mut values, 2, number_of_steps as u64);
    }
}

pub mod general {
    pub struct GeneralStruct {
        array: [u8; 64], // an array o
        slice: [u8; 64], // an array o
        string_entry: str,
    }

    impl GeneralStruct {
        fn print_entry(&self) {
            for entry in &self.slice {
                println!("the entry in the slice is {}", entry);
            }
        }
    }

    fn print_loop(value: &Vec<i32>) {
        let value_clone = value.clone(); // get the cloned value
        for index in &value_clone {
            println!("{} \n", index) // for each value referenced in the index, print out the value index
        }
    }

    fn print_string(s: String) {
        println!("print_String: {}", s);
    }

    fn print_str(s: &str) {
        println!("print_str: {}", s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // lennard-jones double loop test
    #[test]
    fn test_double_loop() {
        let lj_params = lennard_jones_simulations::LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 2,
        };
        // call the double_loop function

        // assert that the result is as expected
        //assert_eq!(result, expected_result)
    }

    #[test]
    fn test_lennard_jones() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let lj_params_new = lennard_jones_simulations::LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 10,
        };
    }

    #[test]
    fn berenden_pull_towards_target() {
        /* mock velocities - T = 300K
           let mut t = 300.0;
           let t0 = 350.0;
           let dt = 0.001;
           let tau = 0.1;
        */
        let t0 = 300.0;

        // Define the new simulation for nve
        let mut new_simulation_md =
            match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
                10, 300.0, 30.0, 10.0, 10.0,
            ) {
                // How to handle errors - we are returning a result or a string
                Ok(atoms) => atoms,
                Err(e) => {
                    eprintln!("Failed to create atoms: {}", e); // Log the error
                    return; // Exit early or handle the error as needed
                }
            };

        lennard_jones_simulations::run_md_nve(&mut new_simulation_md, 1000, 0.5, 10.0, "berendsen");
        let dof = 3 * new_simulation_md.len().saturating_sub(3);
        // compute the final temperature of the system
        let t = lennard_jones_simulations::compute_temperature(&mut new_simulation_md, dof);
        println!("Temperature is {}, and target is {}", t, t0);
        assert!((t - t0).abs() < 5.0, "Temperature should approach target");
    }
}
