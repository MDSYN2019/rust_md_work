
// -- lennard jones potential and force 
pub fn lennard_jones_potential(r: f64, sigma: f64, eps: f64) -> f64 {
    /*
    Return the standard lennard jones function
     */
    if r < 1e-9 { return 0.0; } // Avoid singularity
    let u_ij = 4. * eps * (f64::powi(sigma / r, 12) - f64::powi(sigma / r, 6));
    u_ij
}

pub fn lennard_jones_force(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r < 1e-9 { return 0.0; } // Prevent singularity
    let sr6 = (sigma / r).powi(6);
    24.0 * epsilon * (2.0 * sr6 * sr6 - sr6) / r
}


// -- hard sphere potential and force

pub fn hard_sphere_potential(r: f64, sigma: f64) -> f64 {
    /*    
    Return the hard-sphere potential
     */
    let mut u_ij = 0.0;
    if r < 1e-9 { return 0.0; } // Avoid singularity
    if r < sigma {
        u_ij = 1000000000000000000000.0; // meant to simulate infinity..
    } else {
        u_ij = 0.0; // need to be a floating point number
    }
    u_ij
}
