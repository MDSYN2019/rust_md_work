use sang_md::cell_subdivision::SimulationBox;
use sang_md::lennard_jones_simulations;
use sang_md::lennard_jones_simulations::{
    compute_forces_particles, pbc_update, InitOutput, Particle,
};
use sang_md::molecule::io::{write_gro, write_xtc};

fn snapshot(particles: &[Particle]) -> Vec<Particle> {
    particles.to_vec()
}

fn main() -> Result<(), String> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // A simple point-particle "box of water" style LJ fluid
    let n_particles = 256;
    let temperature = 300.0;
    let mass = 18.015;
    let dt = 0.002;
    let nsteps = 200;
    let box_length = 60.0;

    let mut state = lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
        n_particles,
        temperature,
        mass,
        5.0,
        box_length,
        false,
    )?;

    let particles = match &mut state {
        InitOutput::Particles(particles) => particles,
        InitOutput::Systems(_) => return Err("expected particle initialization".to_string()),
    };

    let simulation_box = SimulationBox {
        x_dimension: box_length,
        y_dimension: box_length,
        z_dimension: box_length,
    };

    let mut subcells = simulation_box.create_subcells(10);
    simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
    compute_forces_particles(particles, box_length, &mut subcells);

    let mut frames: Vec<Vec<Particle>> = Vec::with_capacity(nsteps as usize + 1);
    frames.push(snapshot(particles));

    for step in 0..nsteps {
        let a_old: Vec<_> = particles.iter().map(|p| p.force / p.mass).collect();

        for (p, a) in particles.iter_mut().zip(a_old.iter()) {
            p.velocity += 0.5 * a * dt;
        }

        for p in particles.iter_mut() {
            p.position += p.velocity * dt;
        }

        pbc_update(particles, box_length);

        simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
        compute_forces_particles(particles, box_length, &mut subcells);

        for p in particles.iter_mut() {
            let a_new = p.force / p.mass;
            p.velocity += 0.5 * a_new * dt;
        }

        lennard_jones_simulations::apply_thermostat_berendsen_particles(
            particles,
            temperature,
            0.1,
            dt,
        );

        if step % 10 == 0 {
            frames.push(snapshot(particles));
        }
    }

    write_gro(
        "water_box.gro",
        particles,
        nalgebra::Vector3::new(box_length, box_length, box_length),
        "Water-like point-particle box",
    )?;

    write_xtc(
        "water_box.xtc",
        &frames,
        nalgebra::Vector3::new(box_length, box_length, box_length),
        dt as f32,
    )?;

    println!(
        "Wrote water_box.gro and water_box.xtc for {} particles and {} frames",
        particles.len(),
        frames.len()
    );

    Ok(())
}
