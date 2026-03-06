#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::lennard_jones_force_scalar;

/// High-level Python entry point for building custom MD workflows.
#[pyclass]
pub struct PyMdEngine {
    #[pyo3(get, set)]
    pub sigma: f64,
    #[pyo3(get, set)]
    pub epsilon: f64,
}

#[pymethods]
impl PyMdEngine {
    #[new]
    pub fn new(sigma: f64, epsilon: f64) -> Self {
        Self { sigma, epsilon }
    }

    /// Compute Lennard-Jones force magnitude for a given inter-particle distance.
    pub fn force_at_distance(&self, r: f64) -> f64 {
        lennard_jones_force_scalar(r, self.sigma, self.epsilon)
    }
}

#[pyfunction]
fn lj_force_scalar(r: f64, sigma: f64, epsilon: f64) -> f64 {
    lennard_jones_force_scalar(r, sigma, epsilon)
}

#[pyfunction]
fn python_api_version() -> &'static str {
    "0.1"
}

#[pymodule]
fn sang_md_py(_py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMdEngine>()?;
    module.add_function(wrap_pyfunction!(lj_force_scalar, module)?)?;
    module.add_function(wrap_pyfunction!(python_api_version, module)?)?;
    Ok(())
}
