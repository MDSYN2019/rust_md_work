// https://medium.com/analytics-vidhya/practical-introduction-to-hartree-fock-448fc64c107b

pub mod molecular_hf {
    /*
    https://nznano.blogspot.com/2018/03/simple-quantum-chemistry-hartree-fock.html
     */

    use core::mem::swap;
    use core::mem::take;
    use cute::c; // https://crates.io/crates/cute
    use itertools_num::linspace;
    use num::complex::Complex;
    use std::fs; // filesystems?
    use std::fs::File;

    pub struct HFDataset {
        pub list: Vec<i32>, // this data is private now
        pub average: f32,
        pub name: String,
        pub CoordinatesX: Vec<f32>, // Store X coordinates
        pub CoordinatesY: Vec<f32>, // Store Y coordinates
        pub CoordinatesZ: Vec<f32>, // Store Z Coordinates
                                    //pub filename: String,
    }

    // implement methods for the struct of the HFDataset
    impl HFDataset {
        pub fn OpenStructureFile(filename: &str) {
            // We wish to be able to read the coordinates of the coordinates
            let coordinateInformation = fs::read_to_string(filename);
        }
        // https://nznano.blogspot.com/2018/03/simple-quantum-chemistry-hartree-fock.html
        pub fn psi_STO(minimum: f32, maximum: f32, num: i32) {
            // https://stackoverflow.com/questions/45282970/does-rust-have-an-equivalent-to-pythons-list-comprehension-syntax
            //let mut LinspaceData = c![x.abs(), for x in minimum..maximum];
            let mut LinspaceData = linspace::<f32>(minimum, maximum, num.try_into().unwrap());
            //let mut LinspaceData = (minimum..maximum).filter(|x| x.abs()).collect::<Vec<u32>>();
            let zeta: f64 = 1.0;
            let PI: f64 = 3.14159265358979323846264338327950288;
            let r: f64 = 0.0;
            // Need to convert the values from signed to absolute values in the linspace
            // Rust list comphension equivalents - https://stackoverflow.com/questions/45282970/does-rust-have-an-equivalent-to-pythons-list-comprehension-syntax
            // rename these variables
            let v1 = (0u32..9)
                .filter(|x| x % 2 == 0)
                .map(|x| x.pow(2))
                .collect::<Vec<_>>();

            let v2 = (1..10).filter(|x| x % 2 == 0).collect::<Vec<u32>>();
            let psi_STO = (zeta.powf(3.0) / PI).powf(0.5) * (-1.0 * zeta.powf(r));
        }

        pub fn matchevenodd(&mut self) -> Vec<String> {
            let mut vecString: Vec<String> = Vec::new();
            let buf = self.list.clone();

            for item in buf {
                // Depending on whether the value
                // is even or odd, we will push a different categorical
                // string

                if item % 2 == 0 {
                    vecString.push("Even".to_string());
                } else {
                    vecString.push("Odd".to_string());
                }
            }

            return vecString;
        }

        /*
        Find the largest element within the vector that is in the struct
        */

        pub fn largest_i32(&mut self) -> i32 {
            // https://stackoverflow.com/questions/63353762/cannot-move-out-of-which-is-behind-a-mutable-reference<
            let buf = self.list.clone();
            //swap(&mut buf, &self.list);
            let mut largest = buf[0]; // Take the first item in the mutable list
            for item in buf {
                // loop over the list
                if item > largest {
                    largest = item; // if the item is larger than the previously allocated largest value, then we allocate that value as the largest value
                }
            }
            return largest; // return largest
        }

        pub fn add(&mut self, value: i32) {
            self.list.push(value);
            self.update_average();
        }

        pub fn average(&self) -> f32 {
            self.average
        }

        fn update_average(&mut self) {
            let total: i32 = self.list.iter().sum();
            self.average = total as f32 / self.list.len() as f32;
        }
    }
}
