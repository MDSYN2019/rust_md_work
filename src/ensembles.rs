use rand::Rng;
use rand_dir::{Distribution, Normal};
pub mod ensembles {
    #[derive(Clone, Debug)]
    pub struct ThermostatOption {
        pub target_temperature: f64,
        pub relaxation_time: f64,
    }

    #[derive(Clone, Debug)]
    pub struct BarostatOptions {
        pub target_pressure: f64,
    }

    #[derive(Clone, Debug)]
    pub enum Ensemble {
        Nve,
        Nvt(ThermostatOptions),
        Npt(BarostatOptions),
    }

    impl Ensemble {}
}
