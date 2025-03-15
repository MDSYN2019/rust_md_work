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
    pub struct SimulationBox {
        pub x_dimension: f64,
        pub y_dimension: f64,
        pub z_dimension: f64,
    }

    pub struct MolecularCoordinates {}
}

pub mod lennard_jones_simulations {
    
    use crate::lj_parameters::{lennard_jones_potential, lennard_jones_force, hard_sphere_potential};
    use log::{debug, error, info, trace, warn};
    use nalgebra::{zero, Vector3};
    use ndarray::Array2;
    use rand::prelude::*;
    use rand::Rng;

    #[derive(Clone)]
    pub struct LJParameters {
	// simulation parameters
        pub nsteps: i64, 
        pub n: i64, 
        pub i: i64, 
        pub j: i64,
	// lennard jones parameters
        pub epsilon: f64,
        pub sigma: f64,
	pub epslj: f64,
	// computed interaction values
        pub pot: f64,
        pub rij_sq: f64,
        pub sr2: f64,
        pub sr6: f64,
        pub sr12: f64,
        pub na: i64,
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
	    
            
	    A temperature bath can be acieved by periodically resetting all velocities from the Maxwell-Boltzmann distribution
	    at the desired temperature.
   
            Basically, a function to initialise the velocities of the particles within
            and to try to create a initial velocity pool that is consisten with the given temperature

             */

            let mut rng = rand::thread_rng();
            let sigma_mb = (temp / mass).sqrt();
            let mut velocities: Vec<usize> = Vec::with_capacity(self.lj_parameters.na as usize); // The number of atoms is taken as the capacity for the velocities for certain

            // Compute the random velocities
            let v_x = rng.gen::<f64>() * 2.0 * v_max - v_max;
            let v_y = rng.gen::<f64>() * 2.0 * v_max - v_max;
            let v_z = rng.gen::<f64>() * 2.0 * v_max - v_max;

            println!("velocities are {:?} {:?} {:?}", v_x, v_y, v_z);

            let prob = (-0.5 * (v_x * v_x + v_y * v_y + v_z * v_z)
                / (self.lj_parameters.sigma * self.lj_parameters.sigma))
                .exp();
            let rand_val: f64 = rng.gen();

            // Assign the velocity accroding to the maxwell boltzmann velocity distribution
            if rand_val < prob {
                self.velocity[0] = rng.gen::<f64>() * 2.0 * v_max - v_max;
                self.velocity[1] = rng.gen::<f64>() * 2.0 * v_max - v_max;
                self.velocity[2] = rng.gen::<f64>() * 2.0 * v_max - v_max;
            }
        }

	/*
	The velocity algorithm is a numerical method used in molecular dynamics to integrate
	Newton's equations of motion. It improves upon the standard Velocity algorithm by explicitly
	updating both positions and velocities, making it more accurate for systems where
	velocity dependent forces are dependent
	 */
        fn update_position_verlet(&mut self, acceleration: Vector3<f64>, dt: f64) -> () {
            /*
            Verlet scheme to change the position
            Use the verlet scheme to change the velocity
            */
	    // r = r + dt *  v + 0.5 * dt_sq * a
            self.position[0] += self.velocity[0] * dt + 0.5 * acceleration[0] * dt * dt;
            self.position[1] += self.velocity[1] * dt + 0.5 * acceleration[1] * dt * dt;
            self.position[2] += self.velocity[2] * dt + 0.5 * acceleration[2] * dt * dt;
        }

        fn update_velocity_verlet(
            &mut self,
            acceleration: Vector3<f64>,
            dt: f64,
        ) {
            self.velocity[0] += 0.5 * acceleration[0] * dt;
            self.velocity[1] += 0.5 * acceleration[1] * dt;
            self.velocity[2] += 0.5 * acceleration[2] * dt;
        }

	
	//fn update_position_leapfrog(&mut self, acceleration: Vector3<f64>, dt: f64) -> ()
	//    Another alternative is the so-called half-step 'leapfrog' scheme. The origin of
	//    the name becomes apparent when we write the algorithm
	//    
	//    v(t + 0.5 dt)  = v(t-0.5dt) + dta(t)
	//    r(t + dt) = r(t) + dtv(t + 0.5dt)
	//
	//    The stored quantities are the current positions r(t) and accelerations
	//    a(t) together with the mid-step velocities v(t - 0.5dt)
	//    
	// */
	//{
	//    let mut neg_velocity: Vector3<f64>; 
	//    let mut pos_velocity: Vector3<f64>; 
	//
	//    // get the mod position velocities first 
	//    neg_velocity[0] = self.velocity[0] -  0.5 * acceleration[0] * dt;
	//    neg_velocity[1] = self.velocity[0] -  0.5 * acceleration[0] * dt;
	//    neg_velocity[2] = self.velocity[0] -  0.5 * acceleration[0] * dt;
	//    
	//    pos_velocity[0] = self.velocity[0] +  0.5 * acceleration[0] * dt;
	//    pos_velocity[1] = self.velocity[0] +  0.5 * acceleration[0] * dt;
	//    pos_velocity[2] = self.velocity[0] +  0.5 * acceleration[0] * dt;
	//    
        //    self.position[0] += self.velocity[0] * dt * 0.5 * pos_velocity[0];
        //    self.position[1] += self.velocity[1] * dt + 0.5 * pos_velocity[1];
        //    self.position[2] += self.velocity[2] * dt + 0.5 * pos_velocity[2];
	//    
	//    self.velocity[0] = 0.5 * (neg_velocity[0] + pos_velocity[0]);  
	//    self.velocity[1] = 0.5 * (neg_velocity[1] + pos_velocity[1]);  
	//    self.velocity[2] = 0.5 * (neg_velocity[2] + pos_velocity[2]);  	
	//    
	//}
	
    }
    
    // TODO - implement the Liouville operator??
    
    
    pub fn site_site_energy_calculation(particles: &Vec<Particle>) -> f64 {
         /*
         The coordinates r_ia of a site a in molecule i are stored in the elements r(:, i, a)
	 
         For example, if we have two diatomic molecules, then we have r_1a (site a of molecule 1) and
         r_2a (site a of molecule 2). Each molecule is a diatomic molecule (for example, O2). This means that


	 We already have a set of particles with the lennard jones parameters defined and stored within. Using
	 that data, we need to compute the site_site energy
	 
          */
	 let mut total_energy = 0.0;
	 for i in 0..particles.len()   {
	     for j in (i+1)..particles.len() {
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
	     }
	 }

	 total_energy
     }
    
    pub fn create_atoms_with_set_positions_and_velocities(
        number_of_atoms: i64,
        temp: f64,
        mass: f64,
        v_max: f64,
    ) -> Result<Vec<Particle>, String> {
        /*

	Create number of atoms N atoms with the temperature and mass
	with the mass possible velocity for the atoms also defined
	
         */
	
        let mut rng = thread_rng();
        let mut vector_positions: Vec<Particle> = Vec::new();

        // Create the number of atoms in the system with the system as necessary
        for _ in 0..number_of_atoms {
            let mut particle = Particle {
		// create position for atom 
                position: Vector3::new(
                    // generate x y z position values between -10 and 10
                    rng.gen_range(-10.0..10.0),
                    rng.gen_range(-10.0..10.0),
                    rng.gen_range(-10.0..10.0),
                ),
		// create velocity for atom 
                velocity: Vector3::new(
                    // generate velocity values between -1 and 1
                    rng.gen_range(-1.0..1.0),
                    rng.gen_range(-1.0..1.0),
                    rng.gen_range(-1.0..1.0),
                ),
		// Create the LJ parameters for the atom 
                lj_parameters: (LJParameters {
                    n: 10,
                    i: 0,
                    j: 1,
                    epsilon: 1.0,
                    sigma: 1.0,
                    pot: 0.0,
                    rij_sq: 0.0,
                    sr2: 0.0,
                    sr6: 0.0,
                    sr12: 0.0,
                    epslj: 0.0,
                    nsteps: 1000,
                    na: 1,
                }),	
                force: zero(), // initial force on the atom 
                mass: 0.0, // the mass 
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
	    // update the positions
            particle.update_position_verlet(acceleration, dt);
            // update the velocity
            particle.update_velocity_verlet(acceleration, dt);

            println!(
                "After a iteration step, The position and velocity is {:?} and {:?} ",
                particle.position, particle.velocity
            );
        }
    }

    pub fn compute_forces(mut particles: &mut Vec<Particle>, epsilon: f64, sigma: f64) {
	/*
	 Reset forces before recomputing the force terms
	 */
	for particle in particles.iter_mut() {
	    particle.force = Vector3::zeros();
	}

	let n = particles.len(); // Number of particles in the system

	for i in 0..n {
            for j in (i + 1)..n {
		let r_vec = particles[j].position - particles[i].position;
                let r_ij = (particles[j].position - particles[i].position).norm();
                let force = lennard_jones_force(r_ij, epsilon, sigma); // what are the units?
		let force_vector = (force/ r_ij) * r_vec;
		// Applying forces to both particles (Newton's third law)
                particles[i].force += force_vector; // Apply force to particle i
                particles[j].force -= force_vector; // Apply equal and opposite force to particle j
            }
        }
    }

    pub fn compute_temperature(particles: &mut Vec<Particle>) -> f64 {
        /*
        Compute the current temperature of the system

	To implement a thermostat for a molecular dynamics (MD) simulation in Rust, we need
	to control the temperature of the system by adjusting the velocities of the particles

        This is typically done using methods like Berendson thermostat or velocity rescaling	
         */
        let mut total_kinetic_energy = 0.0;
        let num_particles = particles.len() as f64;

        for particle in particles {
            let velocity_sq = particle.velocity.norm_squared();
            total_kinetic_energy += 0.5 * particle.mass * velocity_sq;
        }
        // T = (2/3) * (KE / N)

	/*
	From the equipartition theorem, the average kinetic energy per particle in a
	system of N particles is given by

	<KE> = 3/2 k_b T

	For a system of N particles, the total kinetic energy is:

	KE = 3/2 * N * k_B * T

	If we assume dimensionless temperature (often used in molecular dynamics simulations), we define:

	T = T / k_B

	Which simplifies the equation:

	T = 2/3 * KE/N

	It defines temperature in terms of kinetic energy. In Molecular Dynamics, temperature is proportional
	to the average kinetic energy per particle
	
	 */
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

    pub fn run_md_nve() {
        /*
        We are now equipt to implement a NVE molecular dynamics simulations.

        define time step and number of steps

         */
        let dt = 0.01;
        let steps = 1000;
        let mut particles = create_atoms_with_set_positions_and_velocities(10, 10.0, 10.0, 10.0);
        // First update the positions
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
    fn test_lennard_jones() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let mut lj_params_new = lennard_jones_simulations::LJParameters {
	    // simulation parameters
	    nsteps: 5,
            n: 3,
            i: 0,
            j: 1,
	    // lennard jones parameters
            epsilon: 1.0,
            sigma: 1.0,
            epslj: 0.0,
	    // computed interaction values
            pot: 0.0,
            rij_sq: 0.0,
            sr2: 0.0,
            sr6: 0.0,
            sr12: 0.0,
            na: 5,
        };

	let mut new_simulation_md =
        match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
            3, 10.0, 10.0, 10.0,
        ) {
            Ok(atoms) => atoms,
            Err(e) => {
                eprintln!("Failed to create atoms: {}", e); // Log the error
                return; // Exit early or handle the error as needed
            }
        };

        // Test case 1 : r = sigma, expected result should be 0
        //assert_eq!(result_1, 0.0);
    }
}
