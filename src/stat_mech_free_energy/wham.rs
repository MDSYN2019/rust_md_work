use super::jarzynski::BOLTZMANN_KCAL_MOL_K;

#[derive(Debug, Clone, PartialEq)]
pub struct UmbrellaWindow {
    pub center: f64,
    pub k_bias: f64,
    pub samples: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhamOptions {
    pub temperature_k: f64,
    pub bin_min: f64,
    pub bin_max: f64,
    pub n_bins: usize,
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for WhamOptions {
    fn default() -> Self {
        Self {
            temperature_k: 300.0,
            bin_min: -5.0,
            bin_max: 5.0,
            n_bins: 200,
            tolerance: 1e-7,
            max_iterations: 10000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct WhamResult {
    pub bin_centers: Vec<f64>,
    pub probabilities: Vec<f64>,
    pub pmf_kcal_mol: Vec<f64>,
    pub free_energy_offsets: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WhamError {
    NoWindows,
    EmptyWindow(usize),
    InvalidBins,
    NonPositiveTemperature,
}

fn harmonic_bias(k_bias: f64, x: f64, center: f64) -> f64 {
    0.5 * k_bias * (x - center).powi(2)
}

pub fn run_wham(windows: &[UmbrellaWindow], options: WhamOptions) -> Result<WhamResult, WhamError> {
    if windows.is_empty() {
        return Err(WhamError::NoWindows);
    }
    if options.temperature_k <= 0.0 {
        return Err(WhamError::NonPositiveTemperature);
    }
    if options.n_bins < 2 || options.bin_max <= options.bin_min {
        return Err(WhamError::InvalidBins);
    }
    for (i, w) in windows.iter().enumerate() {
        if w.samples.is_empty() {
            return Err(WhamError::EmptyWindow(i));
        }
    }

    let beta = 1.0 / (BOLTZMANN_KCAL_MOL_K * options.temperature_k);
    let dz = (options.bin_max - options.bin_min) / options.n_bins as f64;
    let n_windows = windows.len();

    let mut bin_centers = Vec::with_capacity(options.n_bins);
    for b in 0..options.n_bins {
        bin_centers.push(options.bin_min + (b as f64 + 0.5) * dz);
    }

    let mut histogram = vec![vec![0.0_f64; options.n_bins]; n_windows];
    let mut counts_per_window = vec![0.0_f64; n_windows];

    for (i, w) in windows.iter().enumerate() {
        for &x in &w.samples {
            if x >= options.bin_min && x < options.bin_max {
                let idx = ((x - options.bin_min) / dz).floor() as usize;
                if idx < options.n_bins {
                    histogram[i][idx] += 1.0;
                    counts_per_window[i] += 1.0;
                }
            }
        }
    }

    let total_histogram: Vec<f64> = (0..options.n_bins)
        .map(|b| histogram.iter().map(|h| h[b]).sum())
        .collect();

    let mut f_i = vec![0.0_f64; n_windows];
    let mut probabilities = vec![0.0_f64; options.n_bins];

    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..options.max_iterations {
        iterations = iter + 1;

        for b in 0..options.n_bins {
            let x = bin_centers[b];
            let mut denom = 0.0;
            for i in 0..n_windows {
                if counts_per_window[i] > 0.0 {
                    let u = harmonic_bias(windows[i].k_bias, x, windows[i].center);
                    denom += counts_per_window[i] * (beta * (f_i[i] - u)).exp();
                }
            }
            probabilities[b] = if denom > 0.0 {
                total_histogram[b] / denom
            } else {
                0.0
            };
        }

        let norm = probabilities.iter().sum::<f64>() * dz;
        if norm > 0.0 {
            for p in &mut probabilities {
                *p /= norm;
            }
        }

        let mut next_f = vec![0.0_f64; n_windows];
        let mut max_delta: f64 = 0.0;

        for i in 0..n_windows {
            let mut zsum = 0.0;
            for b in 0..options.n_bins {
                let u = harmonic_bias(windows[i].k_bias, bin_centers[b], windows[i].center);
                zsum += probabilities[b] * (-beta * u).exp() * dz;
            }

            if zsum > 0.0 {
                next_f[i] = -(1.0 / beta) * zsum.ln();
            } else {
                next_f[i] = f_i[i];
            }
            max_delta = max_delta.max((next_f[i] - f_i[i]).abs());
        }

        f_i = next_f;
        if max_delta < options.tolerance {
            converged = true;
            break;
        }
    }

    let mut pmf_kcal_mol = vec![f64::INFINITY; options.n_bins];
    for b in 0..options.n_bins {
        if probabilities[b] > 0.0 {
            pmf_kcal_mol[b] = -(1.0 / beta) * probabilities[b].ln();
        }
    }

    if let Some(min_val) = pmf_kcal_mol
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .reduce(f64::min)
    {
        for g in &mut pmf_kcal_mol {
            if g.is_finite() {
                *g -= min_val;
            }
        }
    }

    Ok(WhamResult {
        bin_centers,
        probabilities,
        pmf_kcal_mol,
        free_energy_offsets: f_i,
        iterations,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wham_returns_profile() {
        let windows = vec![
            UmbrellaWindow {
                center: -0.5,
                k_bias: 30.0,
                samples: vec![-0.8, -0.6, -0.55, -0.5, -0.4, -0.35],
            },
            UmbrellaWindow {
                center: 0.0,
                k_bias: 30.0,
                samples: vec![-0.2, -0.1, 0.0, 0.05, 0.1, 0.25],
            },
            UmbrellaWindow {
                center: 0.5,
                k_bias: 30.0,
                samples: vec![0.25, 0.35, 0.45, 0.55, 0.65, 0.8],
            },
        ];

        let opts = WhamOptions {
            temperature_k: 300.0,
            bin_min: -1.0,
            bin_max: 1.0,
            n_bins: 40,
            tolerance: 1e-8,
            max_iterations: 5000,
        };

        let result = run_wham(&windows, opts).expect("WHAM should run");
        assert_eq!(result.bin_centers.len(), 40);
        assert_eq!(result.pmf_kcal_mol.len(), 40);
        assert!(result.probabilities.iter().any(|&p| p > 0.0));
    }
}
