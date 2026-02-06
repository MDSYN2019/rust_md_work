pub mod berendsen {

    pub fn apply_barostat_berendsen_particles(
        particles: &mut Vec<Particle>,
        box_length: &mut f64,
        target_pressure: f64,
        tau_p: f64,
        dt: f64,
        compressability: f64,
    ) -> () {
        // if eithe
        if tau_p <= 0.0 || dt <= 0.0 || compressability <= 0.0 || box_length <= 0.0 {
            return;
        }
        let current_pressure = compute_pressure_particles(particles, *box_length);
    }
}
