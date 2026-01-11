use rand::Rng;
use rand_dir::{Distribution, Normal};

pub mod ensembles {

    #[derive(Clone, Debug)]
    pub struct ThermostatOptions {
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

    impl Ensemble {
        pub fn md_thermostat_andersen() -> () {
            /*

            // basing this from frenkel and smit
                [Initialize system]
                [Compute forces and energy]

                    while t < t_max  do

                    switch = 1
                    integrate-A -- propagates half-step
                    F and E
                switch = 2
                    integrate A -- propagates second-half-step

                    t = t + dt
                    sameple -- sample observables
                     */
        }
    }
}
