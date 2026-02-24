pub struct Jarzynski {
    work_values: Vec<f64>,
}

impl Jarzynski {
    pub fn new() -> Self {
        Jarzynski {
            work_values: Vec::new(),
        }
    }

    pub fn add_work_value(&mut self, work: f64) {
        self.work_values.push(work);
    }

    pub fn calculate_free_energy_difference(&self, beta: f64) -> f64 {
        let sum_exp = self
            .work_values
            .iter()
            .map(|&work| (-beta * work).exp())
            .sum::<f64>();

        let average_exp = sum_exp / self.work_values.len() as f64;

        -1.0 / beta * average_exp.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jarzynski() {
        let mut jarzynski = Jarzynski::new();
        jarzynski.add_work_value(1.0);
        jarzynski.add_work_value(2.0);
        jarzynski.add_work_value(3.0);

        let beta = 1.0; // Inverse temperature (1/kT)
        let delta_f = jarzynski.calculate_free_energy_difference(beta);

        log::info!("Free energy difference: {delta_f:.6}");
    }
}
