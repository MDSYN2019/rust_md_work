/*
Autocorrelation function

The autocorrelation function is one of the most important statistical tools in molecular dynamics, because
it tells you how many
*/

use crate::lennard_jones_simulations::Particle;
use crate::lj_parameters;
use crate::molecule;

#[derive(Debug)]
enum PropertyField {
    Energy,
}

enum PropertyValue {
    Number(f64),
}

fn get_property_value(field: Option<PropertyField>, particle: Particle) -> Option<f64> {
    match field {
        Some(PropertyField::Energy) => Some(particle.energy),
        None => None,
    }
}
pub fn compute_average_val(
    particles: &mut Vec<Particle>,
    field: Option<PropertyField>,
    block_steps: u64,
    number_of_steps: u64,
) -> ()
/*
    Declare the val that we wish to block average; we need to first store the values
    in a dynamic storage (here a vec!) and we then average it through the vec value
    
 */
{
    if number_of_steps % block_steps != 0 {
        eprintln!("We dont have the right number of blocks steps defined for the total number of steps {}", number_of_steps);
    }
    // ensure that the property we want to compute the block average for
    // actually exists within the simulation code
    let mut val: Vec<f64> = vec![];
    for block in particles.chunks(block_steps as usize) {
        println!("Computing property over the block {:?}", block_steps);
        //
    }
}

pub fn autocorrelation_function() -> () {}
