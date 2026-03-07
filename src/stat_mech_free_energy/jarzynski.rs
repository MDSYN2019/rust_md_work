use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub const BOLTZMANN_KCAL_MOL_K: f64 = 0.0019872041;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PullSample {
    pub index: i32,
    pub z: f64,
    pub bilayer_com: f64,
    pub force: f64,
    pub work: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UmbrellaSample {
    pub pull: PullSample,
    pub umbrella_center: f64,
    pub umbrella_k: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeEnergyEstimate {
    pub value: f64,
    pub stdev: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BinResult {
    pub center: i32,
    pub lower: f64,
    pub upper: f64,
    pub raw: FreeEnergyEstimate,
    pub taylor: FreeEnergyEstimate,
    pub alpha: FreeEnergyEstimate,
}

#[derive(Debug)]
pub enum JarzynskiError {
    Io {
        path: String,
        source: std::io::Error,
    },
    InvalidColumns {
        path: String,
        line: usize,
    },
    InvalidNumber {
        path: String,
        line: usize,
    },
    EmptyWorkVector,
    EmptyUmbrellaSamples,
    InvalidTemperature(f64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FreeEnergyMethod {
    Jarzynski,
    UmbrellaSampling,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelectedFreeEnergyEstimate {
    pub method: FreeEnergyMethod,
    pub estimate: FreeEnergyEstimate,
}

pub fn read_pull_file(path: impl AsRef<Path>) -> Result<Vec<PullSample>, JarzynskiError> {
    let path = path.as_ref();
    let path_str = path.display().to_string();
    let file = File::open(path).map_err(|source| JarzynskiError::Io {
        path: path_str.clone(),
        source,
    })?;
    let mut out = Vec::new();
    for (line_idx, line) in BufReader::new(file).lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line.map_err(|source| JarzynskiError::Io {
            path: path_str.clone(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<_> = line.split_whitespace().collect();
        if cols.len() != 5 {
            return Err(JarzynskiError::InvalidColumns {
                path: path_str,
                line: line_no,
            });
        }
        let parse_f64 = |s: &str| s.parse::<f64>().ok();
        let sample = PullSample {
            index: cols[0]
                .parse::<i32>()
                .map_err(|_| JarzynskiError::InvalidNumber {
                    path: path_str.clone(),
                    line: line_no,
                })?,
            z: parse_f64(cols[1]).ok_or_else(|| JarzynskiError::InvalidNumber {
                path: path_str.clone(),
                line: line_no,
            })?,
            bilayer_com: parse_f64(cols[2]).ok_or_else(|| JarzynskiError::InvalidNumber {
                path: path_str.clone(),
                line: line_no,
            })?,
            force: parse_f64(cols[3]).ok_or_else(|| JarzynskiError::InvalidNumber {
                path: path_str.clone(),
                line: line_no,
            })?,
            work: parse_f64(cols[4]).ok_or_else(|| JarzynskiError::InvalidNumber {
                path: path_str.clone(),
                line: line_no,
            })?,
        };
        out.push(sample);
    }
    Ok(out)
}

pub fn read_pull_files(paths: &[impl AsRef<Path>]) -> Result<Vec<PullSample>, JarzynskiError> {
    let mut all = Vec::new();
    for path in paths {
        all.extend(read_pull_file(path)?);
    }
    Ok(all)
}

pub fn raw_jarzynski(
    work: &[f64],
    temperature_k: f64,
) -> Result<FreeEnergyEstimate, JarzynskiError> {
    if work.is_empty() {
        return Err(JarzynskiError::EmptyWorkVector);
    }
    validate_temperature(temperature_k)?;
    let beta = 1.0 / (-BOLTZMANN_KCAL_MOL_K * temperature_k);
    let transformed: Vec<f64> = work.iter().map(|w| (w * beta).exp()).collect();
    let mean = transformed.iter().sum::<f64>() / transformed.len() as f64;
    let value = mean.ln() / beta;
    let scaled: Vec<f64> = transformed.iter().map(|x| x / beta).collect();
    Ok(FreeEnergyEstimate {
        value,
        stdev: stdev_population(&scaled),
    })
}

pub fn taylor_jarzynski(
    work: &[f64],
    temperature_k: f64,
) -> Result<FreeEnergyEstimate, JarzynskiError> {
    if work.is_empty() {
        return Err(JarzynskiError::EmptyWorkVector);
    }
    validate_temperature(temperature_k)?;
    let beta = 1.0 / (-BOLTZMANN_KCAL_MOL_K * temperature_k);
    let mean_w = work.iter().sum::<f64>() / work.len() as f64;
    let mean_w2 = work.iter().map(|w| w * w).sum::<f64>() / work.len() as f64;
    let value = mean_w + (beta / 2.0) * (mean_w2 - mean_w * mean_w);
    Ok(FreeEnergyEstimate {
        value,
        stdev: stdev_population(work),
    })
}

pub fn alpha_jarzynski(
    work: &[f64],
    temperature_k: f64,
) -> Result<FreeEnergyEstimate, JarzynskiError> {
    if work.is_empty() {
        return Err(JarzynskiError::EmptyWorkVector);
    }
    validate_temperature(temperature_k)?;
    let beta = 1.0 / (-BOLTZMANN_KCAL_MOL_K * temperature_k);
    let mean_w = work.iter().sum::<f64>() / work.len() as f64;
    let stdev_w = stdev_population(work);
    let wdiss = 0.5 * beta * stdev_w;
    let alpha = ((15.0 * beta * wdiss).ln()) / ((15.0 * (2.0 * beta * wdiss).exp() - 1.0).ln());
    let bias = wdiss / 10_f64.powf(alpha);

    let transformed: Vec<f64> = work.iter().map(|w| (w * beta).exp()).collect();
    let mean = transformed.iter().sum::<f64>() / transformed.len() as f64;
    let value = mean.ln() / beta - bias;

    let shifted: Vec<f64> = transformed.iter().map(|x| x / beta - bias).collect();
    let _ = mean_w;
    Ok(FreeEnergyEstimate {
        value,
        stdev: stdev_population(&shifted),
    })
}

pub fn umbrella_sampling_free_energy(
    samples: &[UmbrellaSample],
    temperature_k: f64,
) -> Result<FreeEnergyEstimate, JarzynskiError> {
    if samples.is_empty() {
        return Err(JarzynskiError::EmptyUmbrellaSamples);
    }
    validate_temperature(temperature_k)?;
    let beta = 1.0 / (BOLTZMANN_KCAL_MOL_K * temperature_k);

    let unbiased_weights: Vec<f64> = samples
        .iter()
        .map(|sample| {
            let ubias =
                harmonic_bias_energy(sample.pull.z, sample.umbrella_center, sample.umbrella_k);
            (beta * ubias).exp()
        })
        .collect();

    let mean_weight = unbiased_weights.iter().sum::<f64>() / unbiased_weights.len() as f64;
    let value = -mean_weight.ln() / beta;

    Ok(FreeEnergyEstimate {
        value,
        stdev: stdev_population(&unbiased_weights),
    })
}

pub fn choose_free_energy_method(
    method: FreeEnergyMethod,
    work: &[f64],
    umbrella_samples: &[UmbrellaSample],
    temperature_k: f64,
) -> Result<SelectedFreeEnergyEstimate, JarzynskiError> {
    let estimate = match method {
        FreeEnergyMethod::Jarzynski => raw_jarzynski(work, temperature_k)?,
        FreeEnergyMethod::UmbrellaSampling => {
            umbrella_sampling_free_energy(umbrella_samples, temperature_k)?
        }
    };

    Ok(SelectedFreeEnergyEstimate { method, estimate })
}

pub fn run_sample_free_energy_simulation(
    method: FreeEnergyMethod,
    temperature_k: f64,
) -> Result<SelectedFreeEnergyEstimate, JarzynskiError> {
    let sample_work = [0.9, 1.1, 1.0, 1.2];
    let sample_umbrella = vec![
        UmbrellaSample {
            pull: PullSample {
                index: 0,
                z: -0.2,
                bilayer_com: 0.0,
                force: 0.0,
                work: 0.9,
            },
            umbrella_center: -0.1,
            umbrella_k: 5.0,
        },
        UmbrellaSample {
            pull: PullSample {
                index: 1,
                z: 0.0,
                bilayer_com: 0.0,
                force: 0.0,
                work: 1.1,
            },
            umbrella_center: 0.0,
            umbrella_k: 5.0,
        },
    ];

    choose_free_energy_method(method, &sample_work, &sample_umbrella, temperature_k)
}

pub fn compute_bins(
    samples: &[PullSample],
    temperature_k: f64,
) -> Result<Vec<BinResult>, JarzynskiError> {
    if samples.is_empty() {
        return Ok(Vec::new());
    }
    let min_bin = samples
        .iter()
        .map(|s| s.z)
        .fold(f64::INFINITY, f64::min)
        .floor() as i32;
    let max_bin = samples
        .iter()
        .map(|s| s.z)
        .fold(f64::NEG_INFINITY, f64::max)
        .ceil() as i32;

    let mut bins = Vec::new();
    for center in min_bin..=max_bin {
        let lower = center as f64 - 0.5;
        let upper = center as f64 + 0.5;
        let work: Vec<f64> = samples
            .iter()
            .filter(|s| s.z > lower && s.z < upper)
            .map(|s| s.work)
            .collect();
        if work.is_empty() {
            continue;
        }
        bins.push(BinResult {
            center,
            lower,
            upper,
            raw: raw_jarzynski(&work, temperature_k)?,
            taylor: taylor_jarzynski(&work, temperature_k)?,
            alpha: alpha_jarzynski(&work, temperature_k)?,
        });
    }
    Ok(bins)
}

fn stdev_population(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
    var.sqrt()
}

fn validate_temperature(temperature_k: f64) -> Result<(), JarzynskiError> {
    if temperature_k <= 0.0 || !temperature_k.is_finite() {
        return Err(JarzynskiError::InvalidTemperature(temperature_k));
    }
    Ok(())
}

fn harmonic_bias_energy(z: f64, center: f64, k: f64) -> f64 {
    0.5 * k * (z - center).powi(2)
}

impl std::fmt::Display for JarzynskiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            JarzynskiError::Io { path, source } => {
                write!(f, "I/O error while reading `{path}`: {source}")
            }
            JarzynskiError::InvalidColumns { path, line } => write!(
                f,
                "invalid pull row in `{path}` at line {line}: expected 5 columns"
            ),
            JarzynskiError::InvalidNumber { path, line } => {
                write!(f, "invalid numeric value in `{path}` at line {line}")
            }
            JarzynskiError::EmptyWorkVector => write!(f, "no work values were provided"),
            JarzynskiError::EmptyUmbrellaSamples => {
                write!(f, "no umbrella samples were provided")
            }
            JarzynskiError::InvalidTemperature(t) => {
                write!(f, "temperature must be finite and > 0 K, got {t}")
            }
        }
    }
}

impl std::error::Error for JarzynskiError {}
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn parse_pull_file() {
        let mut path = std::env::temp_dir();
        path.push("jarzynski_pull_test.dat");
        fs::write(&path, "0 1.0 0.0 2.0 0.3\n1 1.2 0.0 2.1 0.4\n").unwrap();
        let samples = read_pull_file(&path).expect("parse should succeed");
        assert_eq!(samples.len(), 2);
        let _ = fs::remove_file(path);
    }

    #[test]
    fn computes_estimators() {
        let w = [0.2, 0.3, 0.25, 0.28];
        let raw = raw_jarzynski(&w, 303.0).unwrap();
        let tay = taylor_jarzynski(&w, 303.0).unwrap();
        assert!(raw.value.is_finite());
        assert!(tay.value.is_finite());
        assert!((tay.stdev - stdev_population(&w)).abs() < 1e-12);
    }

    #[test]
    fn computes_umbrella_sampling_free_energy() {
        let samples = vec![
            UmbrellaSample {
                pull: PullSample {
                    index: 0,
                    z: 0.0,
                    bilayer_com: 0.0,
                    force: 0.0,
                    work: 0.9,
                },
                umbrella_center: 0.0,
                umbrella_k: 5.0,
            },
            UmbrellaSample {
                pull: PullSample {
                    index: 1,
                    z: 0.3,
                    bilayer_com: 0.0,
                    force: 0.0,
                    work: 1.1,
                },
                umbrella_center: 0.2,
                umbrella_k: 5.0,
            },
        ];

        let estimate = umbrella_sampling_free_energy(&samples, 300.0).unwrap();
        assert!(estimate.value.is_finite());
        assert!(estimate.stdev.is_finite());
    }

    #[test]
    fn chooser_uses_selected_method() {
        let work = [1.0, 1.2, 1.4];
        let umbrella = vec![UmbrellaSample {
            pull: PullSample {
                index: 0,
                z: 0.0,
                bilayer_com: 0.0,
                force: 0.0,
                work: 1.0,
            },
            umbrella_center: 0.0,
            umbrella_k: 2.0,
        }];

        let jarzynski_choice =
            choose_free_energy_method(FreeEnergyMethod::Jarzynski, &work, &umbrella, 300.0)
                .unwrap();
        assert_eq!(jarzynski_choice.method, FreeEnergyMethod::Jarzynski);

        let umbrella_choice =
            choose_free_energy_method(FreeEnergyMethod::UmbrellaSampling, &work, &umbrella, 300.0)
                .unwrap();
        assert_eq!(umbrella_choice.method, FreeEnergyMethod::UmbrellaSampling);
    }

    #[test]
    fn sample_simulation_runs_for_both_methods() {
        let jarz = run_sample_free_energy_simulation(FreeEnergyMethod::Jarzynski, 300.0).unwrap();
        let umb =
            run_sample_free_energy_simulation(FreeEnergyMethod::UmbrellaSampling, 300.0).unwrap();

        assert!(jarz.estimate.value.is_finite());
        assert!(umb.estimate.value.is_finite());
    }
}
