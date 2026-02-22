/*
Autocorrelation function

The autocorrelation function is one of the most important statistical tools in molecular dynamics, because
it tells you how many
*/


pub fn compute_average_val(
    container_value: &mut Vec<f32>,
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

    for (i, chunk) in container_value.chunks(block_steps as usize).enumerate() {
        let summed_values: f32 = chunk.iter().sum();
        println!(
            "chunk {:?}: {:?} with block sizes as {:?}",
            i,
            summed_values / block_steps as f32,
            block_steps
        );
    }
}

pub fn autocorrelation_function() -> () {}

pub fn radial_distribution_function() -> () {}
