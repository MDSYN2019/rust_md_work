Ipub mod periodic_boundary_conditions {
    /*
    How do we handle periodic boundaries and minimum image convention in a simulation program?

    Define simulation box and write the necessary methods
    to fix the coordiantes when the molecule has periodic boundary issues

     */
    pub struct SimulationBox {
        pub x_dimension: f64,
        pub y_dimension: f64,
        pub z_dimension: f64,
    }

    pub struct MolecularCoordinates {}
}
