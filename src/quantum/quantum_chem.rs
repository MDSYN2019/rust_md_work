// has been separated out
pub mod molecular_polars {
    /*
    Working with the polars library to read in the xyz files from the C++ project of
    quantum chemistry - trying to convert into rust
     */
    use polars::error::PolarsError;
    use polars::prelude::*;
    use std::fs::File;

    pub fn divide(a: f64, b: f64) -> Result<f64, &'static str> {
        if b == 0.0 {
            return Err("Division by zero!");
        }

        Ok(a / b)
    }

    // a struct that can be compared
    //#[derive(Partialeq, PartialOrd)]
    pub struct i_j {
        pub x_i: f64,
        pub y_i: f64,
        pub z_i: f64,
        pub x_j: f64,
        pub y_j: f64,
        pub z_j: f64,
    }

    impl i_j {
        fn rij(&mut self, x_i: f64, x_j: f64, y_i: f64, y_j: f64, z_i: f64, z_j: f64) -> f64 {
            // Possible interatomic distances
            let rij = ((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt();
            rij
        }
    }

    pub fn polars_read_molecular_data_file(readable_file: &str) -> Result<DataFrame, PolarsError> {
        /*
        Here, the standard case is the dataframe, and the err is the polarserror.
        How do we define the dataframe and the polars error here?
         */
        let mut df = CsvReader::from_path(readable_file)?
            .infer_schema(None)
            .has_header(true)
            .finish();
        df
    }

    pub fn polars_read_expression(
        readable_file: &str,
        column: &str,
    ) -> Result<DataFrame, PolarsError> {
        /*
        Select the column we would like from the column entry
         */
        let mut output = polars_read_molecular_data_file(readable_file)
            .expect("Failed to read the molecular data");
        let column_data = output.clone().lazy().select([col(column)]).collect();
        column_data
    }
}

pub mod molecular_structures {
    use core::mem::swap;
    use num::complex::Complex;
    use rand::prelude::*;
    use rand::{random, Rng};
    use rand_distr::{Distribution, Normal};

    /// Linear algebra functionalities
    use ndarray_linalg::norm;
    use std::cmp::Ordering;
    use std::io;
    use std::io::prelude::*;
    use std::io::ErrorKind;
    use std::thread;
    use std::time::Duration;

    // import polars
    use polars::prelude::*;
    use std::fs::File;

    use ndarray::{array, Array1, ArrayView1};

    /*

    The compiler is capable of providing basic implementations of some traits via the #[derive]
    attribute. These traits can still be manually implemented if a more complex behaviour is required.

     */
    pub trait generate_vec {
        fn vectorize(&self) -> Vec<Vec<f32>>;
    }

    #[derive(Debug)] // what are we doing with this debu
    pub struct FILETYPE {
        name: String,
        data: Vec<u8>,
    }

    #[derive(Debug)]
    pub struct gauss_input_coord {
        pub mean: f32,    // mean value
        pub std_dev: f32, // standard deviation
    }

    pub fn polars_read_molecular_data_file(readable_file: &str) -> Result<DataFrame, PolarsError> {
        /*
        A dataframe is a 2-dimensional data structure that is backed by a series, and it could be
        seen as an abstraction of series
        */
        let dataframe = File::open(readable_file.to_string()).expect("could not open file");
        // return type here
        let output = CsvReader::new(dataframe)
            .infer_schema(None)
            .has_header(true)
            .finish();

        output
    }

    fn overlap_integral_readable_file(a: i32, B: i32, RAB2: i32) -> () {} // I think this is the equivalent of a void to be honest

    impl FILETYPE {
        ///
        /// A file interface abel to read a generic molecule mechanics or a force field file for
        /// crunching numbers with this rust program
        ///

        fn new(name: &str) -> FILETYPE {
            /// generate filetype
            /// Generate new FILETYPE with name, but empty vector
            FILETYPE {
                name: String::from(name),
                data: Vec::new(),
            }
        }

        // norms
        pub fn l1_norm(x: ArrayView1<f64>) -> f64 {
            x.fold(0., |acc, elem| acc + elem.abs())
        }

        pub fn l2_norm(x: ArrayView1<f64>) -> f64 {
            x.dot(&x).sqrt()
        }

        //fn normalize(mut x: Array1<f64>) -> Array1<f64> {
        //    let norm = l2_norm(x.view());
        //    x.mapv_inplace(|e| e / norm);
        //   x
        //}

        fn gauss_product(gauss_A: &gauss_input_coord, gauss_B: &gauss_input_coord) -> () {
            let normal_a = Normal::new(gauss_A.mean, gauss_A.std_dev);
            let normal_b = Normal::new(gauss_A.mean, gauss_A.std_dev);
            //let p = normal_a + normal_b;
            //let diff = norm(normal_a - normal_b) ^ 2;
        }

        fn new_with_data(name: &str, data: &Vec<u8>) -> FILETYPE {
            let mut f = FILETYPE::new(name); // generate FILETYPE struct with the string name
            f.data = data.clone();
            f // return f, which is a FILETYPE struct
        }

        #[allow(dead_code)]
        fn read(f: &mut FILETYPE, save_to: &mut Vec<u8>) -> Result<usize, String> {
            /// Result<T, E> -> T is an integer type usize, and E is string. Using String
            let mut tmp = f.data.clone();
            let read_length = tmp.len();
            save_to.reserve(read_length); // Ensures that there is sufficient space to fit the incoming data
            save_to.append(&mut tmp);
            Ok(read_length) // return read length. Otherwise, return a string as we are returning a Result
        }
    }

    // main structure
    pub struct LinesWithEndings<'a> {
        input: &'a str, // ????
    }

    impl<'a> LinesWithEndings<'a> {
        pub fn from(input: &'a str) -> LinesWithEndings<'a> {
            LinesWithEndings { input: input }
        }
    }

    // Implement method
    impl<'a> Iterator for LinesWithEndings<'a> {
        type Item = &'a str;
        #[inline]
        fn next(&mut self) -> Option<&'a str> {
            if self.input.is_empty() {
                return None;
            }
            let split = self
                .input
                .find('\n')
                .map(|i| i + 1)
                .unwrap_or(self.input.len());
            let (line, rest) = self.input.split_at(split);
            self.input = rest;
            Some(line)
        }
    }

    // functions within the mods

    fn one_in(denominator: u32) -> bool {
        thread_rng().gen_ratio(1, denominator)
    }

    fn open(f: FILETYPE) -> Result<FILETYPE, String> {
        if one_in(10_000) {
            let err_msg = String::from("Permission denied");
            return Err(err_msg);
        }
        Ok(f)
    }
    fn close(f: FILETYPE) -> Result<FILETYPE, String> {
        if one_in(100_000) {
            // Once in 10000 executions, return an error
            let err_msg = String::from("Interrupted by signal!");
            return Err(err_msg);
        }
        Ok(f)
    }
    // specific coordinate operation
    pub struct Coordinates {
        pub X: f32, // x coordinates
        pub Y: f32, // y coordinates
    }

    impl generate_vec for Coordinates {
        fn vectorize(&self) -> Vec<Vec<f32>> {
            let mut vector_of_vectors: Vec<Vec<f32>> = Vec::new(); // initialize new vector of vectors
            let mut inner_vector: Vec<f32> = Vec::new(); // initialize new vector
            inner_vector.push(self.X);
            inner_vector.push(self.Y);
            //vector_of_vectors.push(v); // dont have v at the moment - dont know where this is for now ..
            vector_of_vectors
        }
    }
}

pub mod self_consistent_field {
    /*
    implementing the self-consistent field theory implmentation for extracting the wavefunction

    The HF-self_consistent_field is often formalised in such abstractness it does not become clear how
    the method really works

    Need to knows:

    (i) Dirac notation

    (ii) The stationary wave functions of the coulomb central-field problem

    (iii) The variational principle for bound states

    (iv) The idea of completeness

     */
    use super::*;
    use core::mem::{swap, take};
    use cute::c; // https://crates.io/crates/cute
    use itertools_num::linspace;
    use kdam::tqdm; // tqdm - rust version!
    use num::complex::Complex;
    use polars::prelude::*;
    use std::fs; // filesystems?
    use std::fs::File; // import dataframe from here
                       // implementing Santra and Obermeyer

    pub struct scf_vals {
        pub total_energy: f64,
    }

    pub struct atomic_parameters {
        pub atomic_number: f64,
        pub num_iter: i32,
        // define the six coulomb integrals
        pub eps_1: f64,
        pub eps_2: f64,
        pub I_1111: f64,
        pub I_1112: f64,
        pub I_1122: f64,
        pub I_1212: f64,
        pub I_1222: f64,
        pub I_2222: f64,
    }
    // define types

    type scf_val_type = scf_vals;
    type atomic_parameters_type = atomic_parameters;

    // make default values for atomic_parameters
    impl Default for atomic_parameters {
        fn default() -> atomic_parameters {
            // return type atomic parameters
            atomic_parameters {
                // initialized values for the HF parameters
                atomic_number: 2.0,
                num_iter: 100,
                eps_1: 0.0,
                eps_2: 0.0,
                I_1111: 0.0,
                I_1112: 0.0,
                I_1122: 0.0,
                I_1212: 0.0,
                I_1222: 0.0,
                I_2222: 0.0,
            }
        }
    }

    pub struct input_files {
        /// place to store the polar input files
        f_theta: f64,
        f_total_energy: f64,
        f_orbital_energy: f64,
    }

    impl atomic_parameters {
        fn read(readable_file: &str) {
            let dataframe = File::open(readable_file.to_string()).expect("could not open file");
        }

        fn assert_values(&mut self) {
            assert_eq!(self.eps_1, self.eps_1);
            assert_eq!(self.eps_2, self.eps_2);
            assert_eq!(self.I_1111, self.I_1111);
        }

        pub fn compute_I_values(&mut self) {
            // set the values of the six Coulomb integrals
            self.eps_1 = -1. * f64::powi(self.atomic_number, 2) / 2.;
            self.eps_2 = -1. * f64::powi(self.atomic_number, 2) / 8.;
            self.I_1111 = (5. / 8.) * self.atomic_number;
            self.I_1112 =
                f64::powi(2., 12) * (2_f64).sqrt() / 27. / f64::powi(7., 4) * self.atomic_number;
            self.I_1122 = f64::powi(16. / 9., 3) * self.atomic_number;
            self.I_1222 =
                f64::powi(2., 9) * (2_f64).sqrt() / 27. / f64::powi(5., 5) * self.atomic_number;
            self.I_1212 = f64::powi((17 / 3) as f64, 4) * self.atomic_number;
            self.I_2222 = f64::powi((77 / 2) as f64, 9) * self.atomic_number;
        }

        pub fn read_lines(&mut self, filename: &str) -> Vec<String> {
            fs::read_to_string(filename)
                .unwrap()
                .lines()
                .map(String::from)
                .collect()
        }

        fn set_theta_energy_orbital(
            &mut self,
            theta: f64,
            total_energy: f64,
            orbital_energy: f64,
        ) -> input_files {
            /// define another type and return the struct
            /// more lines
            input_files {
                f_theta: theta,
                f_total_energy: total_energy,
                f_orbital_energy: orbital_energy,
            }
        }

        pub fn compute_two_electron_energy(&mut self) -> () {
            /// some commentary here like a python docstring
            // public function
            let mut total_energy = 0.0;
            let mut theta_energy_orbital = self.set_theta_energy_orbital(0.0, 0.0, 0.0); // private function within the mod called
            let mut c_1: f64;
            let mut c_2: f64;
            // fock matrix components
            let mut F_11: f64 = 0.;
            let mut F_12: f64 = 0.;
            let mut F_21: f64 = 0.;
            let mut F_22: f64 = 0.;
            let mut orbital_energy: f64 = 0.;

            // loop over the number of iterations defined
            for entry in tqdm!(0..self.num_iter) {
                c_1 = theta_energy_orbital.f_theta.cos();
                c_2 = theta_energy_orbital.f_theta.sin();
                total_energy = (2.0 as f64)
                    * (f64::powi(self.eps_1 * c_1, 2) + f64::powi(self.eps_2 * c_2, 2))
                    + (f64::powi(self.I_1111 * c_1, 4))
                    + ((4. * f64::powi(self.I_1112 * c_1, 3)) * c_2)
                    + (f64::powi(2. * (2. * self.I_1122 + self.I_1212) * c_1, 2)
                        * f64::powi(2. * c_2, 2))
                    + (4. * f64::powi(self.I_1222 * c_1 * c_2, 3))
                    + f64::powi(self.I_2222 * c_2, 4);

                log::info!(
                    "SCF iter {entry:>4} | E_total={total_energy:.8} c1={c_1:.6} c2={c_2:.6}"
                );

                // Computing the fock matrix 1,1 th element?
                F_11 = f64::powi(self.eps_1 + self.I_1111 * c_1, 2)
                    + 2. * (self.I_1112 * c_1 * c_2)
                    + f64::powi(self.I_1212 * c_2, 2);
                // 1,2 th element?
                F_12 = f64::powi(self.I_1112 * c_1, 2)
                    + (2. * self.I_1212 * c_1 * c_2)
                    + f64::powi(self.I_1222 * c_2, 2);
                // 2,1 th element?
                F_21 = F_12;
                F_22 = self.eps_2
                    + f64::powi(self.I_1212 * c_1, 2)
                    + (2. * self.I_1222 * c_1 * c_2)
                    + f64::powi(self.I_2222 * c_2, 2);

                // calculation of the lower of the two roots of the characteristic
                // polynomial of the Fock matrix, using the quadratic formula
                orbital_energy =
                    0.5 * (F_11 + F_22) - (0.25 * f64::powi(F_11 - F_22, 2) + F_12 * F_21).sqrt();
                // arctan
                theta_energy_orbital.f_theta = ((orbital_energy - F_11) / F_12).atan();
            }
        }
    }
}
