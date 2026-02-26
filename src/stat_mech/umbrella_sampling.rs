#[derive(Debug, Clone)]
pub struct UmbrellaWindow {
    pub center: f64,
    pub force_constant: f64,
    samples: Vec<f64>,
}

impl UmbrellaWindow {
    pub fn new(center: f64, force_constant: f64) -> Self {
        Self {
            center,
            force_constant,
            samples: Vec::new(),
        }
    }

    pub fn add_sample(&mut self, reaction_coordinate: f64) {
        self.samples.push(reaction_coordinate);
    }

    pub fn add_samples<I>(&mut self, coordinates: I)
    where
        I: IntoIterator<Item = f64>,
    {
        self.samples.extend(coordinates);
    }

    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    fn bias_potential(&self, reaction_coordinate: f64) -> f64 {
        let displacement = reaction_coordinate - self.center;
        0.5 * self.force_constant * displacement * displacement
    }
}

#[derive(Debug, Clone)]
pub struct UmbrellaSampling {
    windows: Vec<UmbrellaWindow>,
    min_coordinate: f64,
    max_coordinate: f64,
    n_bins: usize,
    temperature: f64,
    boltzmann_constant: f64,
}

#[derive(Debug, Clone)]
pub struct PmfProfile {
    pub bin_centers: Vec<f64>,
    pub free_energies: Vec<f64>,
}

impl UmbrellaSampling {
    pub fn new(
        windows: Vec<UmbrellaWindow>,
        min_coordinate: f64,
        max_coordinate: f64,
        n_bins: usize,
        temperature: f64,
    ) -> Self {
        Self {
            windows,
            min_coordinate,
            max_coordinate,
            n_bins,
            temperature,
            boltzmann_constant: 1.0,
        }
    }

    pub fn with_boltzmann_constant(mut self, boltzmann_constant: f64) -> Self {
        self.boltzmann_constant = boltzmann_constant;
        self
    }

    pub fn calculate_pmf(&self) -> Result<PmfProfile, String> {
        if self.windows.is_empty() {
            return Err("Umbrella sampling requires at least one window".to_string());
        }
        if self.n_bins < 2 {
            return Err("Umbrella sampling requires at least 2 bins".to_string());
        }
        if self.max_coordinate <= self.min_coordinate {
            return Err("max_coordinate must be larger than min_coordinate".to_string());
        }
        if self.temperature <= 0.0 {
            return Err("temperature must be larger than zero".to_string());
        }

        let beta = 1.0 / (self.boltzmann_constant * self.temperature);
        let bin_width = (self.max_coordinate - self.min_coordinate) / self.n_bins as f64;

        let mut total_samples = 0usize;
        let mut weighted_probability = vec![0.0; self.n_bins];

        for window in &self.windows {
            if window.samples.is_empty() {
                continue;
            }

            total_samples += window.samples.len();

            for &coordinate in &window.samples {
                if let Some(bin_index) = self.coordinate_to_bin(coordinate) {
                    let bias_energy = window.bias_potential(coordinate);
                    weighted_probability[bin_index] += (beta * bias_energy).exp();
                }
            }
        }

        if total_samples == 0 {
            return Err("No samples found in umbrella windows".to_string());
        }

        let normalization: f64 = weighted_probability.iter().sum();
        if normalization <= 0.0 {
            return Err("Unable to normalize umbrella histogram".to_string());
        }

        let mut free_energies = Vec::with_capacity(self.n_bins);
        let mut bin_centers = Vec::with_capacity(self.n_bins);

        for (bin_index, &weight) in weighted_probability.iter().enumerate() {
            let center = self.min_coordinate + (bin_index as f64 + 0.5) * bin_width;
            bin_centers.push(center);

            if weight > 0.0 {
                let probability = weight / normalization;
                free_energies.push(-(1.0 / beta) * probability.ln());
            } else {
                free_energies.push(f64::INFINITY);
            }
        }

        let baseline = free_energies
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .fold(f64::INFINITY, f64::min);

        if baseline.is_finite() {
            for value in &mut free_energies {
                if value.is_finite() {
                    *value -= baseline;
                }
            }
        }

        Ok(PmfProfile {
            bin_centers,
            free_energies,
        })
    }

    fn coordinate_to_bin(&self, coordinate: f64) -> Option<usize> {
        if coordinate < self.min_coordinate || coordinate >= self.max_coordinate {
            return None;
        }

        let fraction =
            (coordinate - self.min_coordinate) / (self.max_coordinate - self.min_coordinate);
        let mut index = (fraction * self.n_bins as f64).floor() as usize;
        if index >= self.n_bins {
            index = self.n_bins - 1;
        }
        Some(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn umbrella_reweighting_recovers_flat_profile() {
        let mut left_window = UmbrellaWindow::new(-1.0, 0.0);
        let mut right_window = UmbrellaWindow::new(1.0, 0.0);

        for z_index in 0..=40 {
            let z = -2.0 + z_index as f64 * 0.1;

            for _ in 0..20 {
                left_window.add_sample(z);
                right_window.add_sample(z);
            }
        }

        let sampler = UmbrellaSampling::new(vec![left_window, right_window], -2.0, 2.0, 40, 1.0)
            .with_boltzmann_constant(1.0);

        let profile = sampler
            .calculate_pmf()
            .expect("PMF should be computed for valid windows");

        let finite_values: Vec<f64> = profile
            .free_energies
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .collect();

        let mean = finite_values.iter().sum::<f64>() / finite_values.len() as f64;
        let variance = finite_values
            .iter()
            .map(|value| {
                let delta = value - mean;
                delta * delta
            })
            .sum::<f64>()
            / finite_values.len() as f64;

        assert!(variance.sqrt() < 3.0);
    }

    #[test]
    fn rejects_empty_window_data() {
        let sampler =
            UmbrellaSampling::new(vec![UmbrellaWindow::new(0.0, 5.0)], -1.0, 1.0, 20, 1.0);

        let err = sampler
            .calculate_pmf()
            .expect_err("sampler should fail without samples");

        assert!(err.contains("No samples"));
    }
}
