/*
Code based on the following paper:

https://apt.scitation.org/doi/abs/10.1119/10.0002644?journalCode=ajp

Title: A first encounter with the Hartree-fock self-consistent field method

*/

// Declaring the modules and importing from it
//pub mod stat_mech;

extern crate assert_type_eq;
mod lj_parameters;

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
    use nalgebra::{zero, Vector3};

    pub struct SimulationBox {
        pub x_dimension: f64,
        pub y_dimension: f64,
        pub z_dimension: f64,
    }

    impl SimulationBox {
        fn cell_subdivison(&self, n_cells: i64) -> () {
            /*
            cell subdivision provides a mean for organizing the information about atom positions
            into a form that avoids most of the unnecessary work and reduces the computational effort to O(N_m) level.

            linked lists are used to associate atoms with the cells in which they reside at any given instant. A separate list is required for each cell.

             */

            let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);
            let cell_size = box_size / (n_cells as f64);
        }
    }
    pub struct MolecularCoordinates {}
}

pub mod lennard_jones_simulations {

    use crate::lj_parameters::{
        hard_sphere_potential, lennard_jones_force, lennard_jones_potential,
    };
    use log::{debug, error, info, trace, warn};
    use nalgebra::{zero, Vector3};
    use ndarray::Array2;
    use rand::prelude::*;
    use rand::Rng;

    #[derive(Clone)]
    pub struct LJParameters {
        // lennard jones parameters and the number of atoms that we have of that parameter
        pub epsilon: f64,
        pub sigma: f64,
        pub number_of_atoms: i32,
    }

    #[derive(Clone)]
    pub struct Particle {
        position: Vector3<f64>,
        velocity: Vector3<f64>,
        force: Vector3<f64>,
        lj_parameters: LJParameters,
        mass: f64,
    }

    impl Particle {
        fn distance(&self, other: &Particle) -> f64 {
            // Compute the distance between two particles
            (self.position - other.position).norm()
        }

        fn maxwellboltzmannvelocity(&mut self, temp: f64, mass: f64, v_max: f64) -> () {
            /*
                ---------------------
                More information here:
                ---------------------
                https://scicomp.stackexchange.com/questions/19969/how-do-i-generate-maxwell-boltzmann-variates-using-a-uniform-distribution-random

                A temperature bath can be achieved by periodically resetting all velocities from the Maxwell-Boltzmann distribution
                at the desired temperature

                A function to initialize the velocities of the particles within and to try to create an initial
                velocity pool that is consistent with the given temperature

            ---

            loop over the total velocities

                 */
            let mut rng = rand::rng();
            let sigma_mb = (temp / mass).sqrt();
            let mut velocities: Vec<usize> =
                Vec::with_capacity(self.lj_parameters.number_of_atoms as usize); // The number of atoms is taken as the capacity for the velocities for certain

            // Compute the random velocities -
            let v_x = rng.random::<f64>() * 2.0 * v_max - v_max; // velocity x component
            let v_y = rng.random::<f64>() * 2.0 * v_max - v_max; // velocity y component
            let v_z = rng.random::<f64>() * 2.0 * v_max - v_max; // velocity z componeont

            println!("velocities are {:?} {:?} {:?}", v_x, v_y, v_z);

            let prob = (-0.5 * (v_x * v_x + v_y * v_y + v_z * v_z)
                / (self.lj_parameters.sigma * self.lj_parameters.sigma))
                .exp(); // Not sure what this is doing - getting the probability distribution?

            let rand_val: f64 = rng.random();

            // Assign the velocity according to the maxwell boltzmann velocity distribution
            if rand_val < prob {
                self.velocity[0] = rng.random::<f64>() * 2.0 * v_max - v_max;
                self.velocity[1] = rng.random::<f64>() * 2.0 * v_max - v_max;
                self.velocity[2] = rng.random::<f64>() * 2.0 * v_max - v_max;
            }
        }

        fn update_position_verlet(&mut self, acceleration: Vector3<f64>, dt: f64) -> () {
            /*
            Verlet scheme to change the position
            Use the verlet scheme to change the velocity
            */
            self.position[0] += self.velocity[0] * dt + 0.5 * acceleration[0] * dt * dt;
            self.position[1] += self.velocity[1] * dt + 0.5 * acceleration[1] * dt * dt;
            self.position[2] += self.velocity[2] * dt + 0.5 * acceleration[2] * dt * dt;
        }

        fn update_velocity_verlet(
            &mut self,
            old_acceleration: Vector3<f64>,
            new_acceleration: Vector3<f64>,
            dt: f64,
        ) {
            self.velocity[0] += 0.5 * (old_acceleration[0] + new_acceleration[0]) * dt;
            self.velocity[1] += 0.5 * (old_acceleration[1] + new_acceleration[1]) * dt;
            self.velocity[2] += 0.5 * (old_acceleration[2] + new_acceleration[2]) * dt;
        }
    }

    pub fn site_site_energy_calculation(particles: &Vec<Particle>) -> f64 {
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

                let sigma_i = particles[i].lj_parameters.sigma;
                let epsilon_i = particles[i].lj_parameters.epsilon;
                let sigma_j = particles[j].lj_parameters.sigma;
                let epsilon_j = particles[j].lj_parameters.epsilon;
                // Using Lorentz-Bethelot mixing rules
                let computed_sigma = (sigma_i + sigma_j) / 2.0;
                let computed_epsilon = (epsilon_i + epsilon_j).sqrt();
                let r_vec = particles[j].position - particles[i].position;
                let r = r_vec.norm();
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
    ) -> Result<Vec<Particle>, String> {
        /*

            Create N atoms with temperature and mass

            Here, we are going with the assumption that we are creating a simulation box
            that is cubic, meaning that we will get a (min + max) * (min + max) *  (min+max) volume
        system for all the molecules

             */
        let mut vector_positions: Vec<Particle> = Vec::new();
        let mut rng = rand::rng();
        // Create the number of atoms in the system with the system as necessary
        for _ in 0..number_of_atoms {
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
                mass: 0.0,     // the mass
            };

            // Reset the positions to the maxwell boltzmann distibution of velocities
            particle.maxwellboltzmannvelocity(temp, mass, v_max);
            // push those values into the vector
            vector_positions.push(particle); // push the newly assigned particle into the positions
        }
        Ok(vector_positions)
    }

    pub fn run_verlet_update(
        particles: &mut Vec<Particle>,
        acceleration: Vector3<f64>,
        dt: f64,
    ) -> () {
        /*
        Update the position and velocity of the particle using the verlet scheme
         */

        // update the position
        for particle in particles.iter_mut() {
            println!(
                "The original position and velocity is {:?} and {:?} ",
                particle.position, particle.velocity
            );
            particle.update_position_verlet(acceleration, dt);
            // update the velocity
            particle.update_velocity_verlet(acceleration, acceleration, dt);

            println!(
                "After a iteration step, the position and velocity is {:?} and {:?} ",
                particle.position, particle.velocity
            );
        }
    }

    pub fn compute_forces(mut particles: &mut Vec<Particle>, epsilon: f64, sigma: f64) {
        // TODO
        let n = particles.len(); // number of particles in the system

        for i in 0..n {
            for j in (i + 1)..n {
                let r_ij = (particles[j].position - particles[i].position).norm();
                let force = lennard_jones_potential(r_ij, epsilon, sigma);
                //particles[i].force += force; // Apply force to particle i
                //particles[j].force -= force; // Apply equal and opposite force to particle j
            }
        }
    }

    pub fn compute_temperature(particles: &mut Vec<Particle>) -> f64 {
        /*
        Compute the current temperature of the system

        To implement a thermostat for a molecular dynamics (MD) simulation in Rust, we need
        to control the temperature of the system by adjusting the velocities of the particles

        This is typically done using methods like Berendson thermostat or velocity rescaling

        ---

        T = 2/3 * (KE/N) - this is the fundamental relation in classical statistical mechanics that connects
        the temperature of a system with its kinetic energy per particle

         */
        let mut total_kinetic_energy = 0.0;
        let num_particles = particles.len() as f64;

        for particle in particles {
            let velocity_sq = particle.velocity.norm_squared();
            total_kinetic_energy += 0.5 * particle.mass * velocity_sq;
        }
        // T = (2/3) * (KE / N)
        (2.0 / 3.0) * (total_kinetic_energy / num_particles)
    }

    pub fn apply_thermostat(particles: &mut Vec<Particle>, target_temperature: f64) {
        let current_temperature = compute_temperature(particles);
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

    pub fn pbc_update(particles: &mut Vec<Particle>, box_length: f64) -> Vec<Particle> {
        // Automatically correct the particle position based on
        // the periodic boundary conditions
        for i in 0..particles.len() {
            particles[i].position[0] = particles[i].position[0] % box_length;
            particles[i].position[1] = particles[i].position[1] % box_length;
            particles[i].position[2] = particles[i].position[2] % box_length;
        }

        particles.clone()
    }

    pub fn run_md_nve(number_of_steps: i32, dt: f64) {
        /*
        We are now equipt to implement a NVE molecular dynamics simulations.

        define time step and number of steps

         */

        let lj_params_new = LJParameters {
            epsilon: 1.0,
            sigma: 4.0,
            number_of_atoms: 2,
        };

        let mut new_simulation_md =
            match create_atoms_with_set_positions_and_velocities(3, 10.0, 10.0, 10.0, 10.0) {
                // How to handle errors - we are returning a result or a string
                Ok(atoms) => atoms,
                Err(e) => {
                    eprintln!("Failed to create atoms: {}", e); // Log the error
                    return; // Exit early or handle the error as needed
                }
            };

        // loop over the total system for number_of_steps
        for i in 0..number_of_steps {
            let mut updated_sim = pbc_update(&mut new_simulation_md, 20.0);
            compute_forces(&mut updated_sim, lj_params_new.epsilon, lj_params_new.sigma);
            // update velocities using the verlet format
            run_verlet_update(&mut updated_sim, Vector3::new(0.01, 0.01, 0.01), 0.05);
        }
    }
}

pub mod general {
    pub struct GeneralStruct {
        // borrwed a slice of an array
        array: [u8; 64], // an array o
        //slice: &array,
        slice: [u8; 64], // an array o
        string_entry: str,
    }

    impl GeneralStruct {
        ///
        ///
        fn print_entry(&self) {
            for entry in &self.slice {
                println!("the entry in the slice is {}", entry);
            }
        }

        //fn split_string(&self) {
        //    for word in &self.string_entry.chars {
        //       println!("{}", word);
        //   }
    }

    fn print_loop(value: &Vec<i32>) {
        let value_clone = value.clone(); // get the cloned value
        for index in &value_clone {
            println!("{} \n", index) // for each value referenced in the index, print out the value index
        }
    }
    // Vec inherits the methods of slices, because we can obtain a slice reference
    // from a vector

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
    /*
    test the self-consistent field theory part here
     */
    #[test]
    fn test_self_consistent_field() {
        //assert_eq!(
        //    self_consistent_field::atomic_parameters::default().atomic_number,
        //    self_consistent_field::atomic_parameters::default().atomic_number
        //);
        //assert_eq! (
        //        self_consistent_field::atomic_parameters::default()
        //        self_consistent_field::atomic_parameters::default()
        //);
    }

    /*
    test implementation of scf_Vals
     */
    //#[test]
    //fn test_energy() {
    //    let scf_entry = self_consistent_field::scf_vals { total_energy: 32.0 };
    //    let scf_entry_wrong = self_consistent_field::scf_vals { total_energy: 20.0 };
    //    //assert_eq!(scf_entry, scf_entry_wrong);
    //}
    //
    //#[test]
    //fn test_atomic_parameters() {
    //    /*
    //    Description of what we are testing here
    //    */
    //    let mut atom_parameters = self_consistent_field::atomic_parameters {
    //        ..Default::default()
    //    };
    //    atom_parameters.compute_I_values();
    //    atom_parameters.print_f_values();
    //    atom_parameters.compute_two_electron_energy();
    //}

    #[test]
    fn test_outer_product() {
        let a = vec![1, 2, 3];
        let b = vec![4, 5, 6];

        let result = tensors::outer_product(&a, &b, 0);
        //assert!(result.is::<Vec<Vec<i32>>());
    }

    #[test]
    fn test_polars() {
        let mut ex = molecular_polars::polars_read_molecular_data_file(
            "/home/sang/Desktop/git/ProgrammingProjects/Project#01/input/benzene.dat",
        );

        let data: Vec<&str> = vec![
            "6    0.000000000000     0.000000000000     0.000000000000",
            "6    0.000000000000     0.000000000000     2.616448463377",
            // Add more strings here
        ];

        //let mut newcontents =
        //    read_lines("/home/sang/Desktop/git/ProgrammingProjects/Project#01/input/benzene.dat");
    }

    // lennard-jones double loop test
    #[test]
    fn test_double_loop() {
        let mut lj_params = LennardJonesSimulations::LJParameters {
            eps: 1.0,
            sigma: 1.0,
            number_of_atoms: 2,
        };
        // call the double_loop function
        let result = lj_params.double_loop();

        // assert that the result is as expected
        //assert_eq!(result, expected_result)
    }

    #[test]
    fn test_lennard_jones() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let mut lj_params_new = LennardJonesSimulations::LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 3,
        };

        // Test case 1 : r = sigma, expected result should be 0
        let r_1 = sigma;
        let result_1 = lj_params_new.lennard_jones(sigma, r_1, epsilon);
        //assert_eq!(result_1, 0.0);
    }
}
