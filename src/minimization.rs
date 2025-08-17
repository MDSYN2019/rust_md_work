/*
Energy minimization in molecular dynamics is a preprocessing step
used to find nearby local minimum of potential energy surface

The idea is to remove steric clashes, highly strained geometries,
or unrealistic bond lengths/angbles before running time-evolution simulations
*/

use argmin::prelude::*;
use argmin::solver::quasinewton::BFGS;
