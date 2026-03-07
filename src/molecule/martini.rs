use crate::lennard_jones_simulations::{LJParameters, Particle};
use crate::molecule::molecule::{Angle, Bond, Dihedral, System};
use nalgebra::Vector3;
use std::collections::HashMap;
use std::fs;

#[derive(Clone, Debug, Default)]
pub struct MartiniAtomType {
    pub name: String,
    pub mass: f64,
    pub charge: f64,
    pub sigma: Option<f64>,
    pub epsilon: Option<f64>,
    pub c6: Option<f64>,
    pub c12: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct MartiniAtom {
    pub index: usize,
    pub type_name: String,
    pub charge: f64,
    pub mass: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct MartiniBond {
    pub atom1: usize,
    pub atom2: usize,
    pub length: f64,
    pub force_constant: f64,
}

#[derive(Clone, Debug)]
pub struct MartiniAngle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub theta0_deg: f64,
    pub force_constant: f64,
}

#[derive(Clone, Debug)]
pub struct MartiniDihedral {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub phase_deg: f64,
    pub force_constant: f64,
    pub multiplicity: usize,
}

#[derive(Clone, Debug, Default)]
pub struct MartiniForceField {
    pub molecule_name: Option<String>,
    pub atom_types: HashMap<String, MartiniAtomType>,
    pub atoms: Vec<MartiniAtom>,
    pub bonds: Vec<MartiniBond>,
    pub angles: Vec<MartiniAngle>,
    pub dihedrals: Vec<MartiniDihedral>,
}

impl MartiniForceField {
    pub fn parse_str(contents: &str) -> Result<Self, String> {
        let mut ff = MartiniForceField::default();
        let mut section = String::new();

        for (line_number, raw_line) in contents.lines().enumerate() {
            let line = raw_line.split(';').next().unwrap_or("").trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.starts_with('[') && line.ends_with(']') {
                section = line[1..line.len() - 1].trim().to_ascii_lowercase();
                continue;
            }

            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            match section.as_str() {
                "moleculetype" => {
                    if ff.molecule_name.is_none() {
                        ff.molecule_name = Some(tokens[0].to_string());
                    }
                }
                "atomtypes" => {
                    let atom_type = parse_atomtype(&tokens)
                        .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                    ff.atom_types.insert(atom_type.name.clone(), atom_type);
                }
                "atoms" => {
                    ff.atoms.push(
                        parse_atom(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "bonds" => {
                    ff.bonds.push(
                        parse_bond(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "angles" => {
                    ff.angles.push(
                        parse_angle(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "dihedrals" => {
                    ff.dihedrals.push(
                        parse_dihedral(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                _ => {}
            }
        }

        if ff.atoms.is_empty() {
            return Err("martini input did not contain an [ atoms ] section".to_string());
        }

        Ok(ff)
    }

    pub fn read_itp(path: &str) -> Result<Self, String> {
        let contents = fs::read_to_string(path)
            .map_err(|e| format!("failed to read martini file at '{path}': {e}"))?;
        Self::parse_str(&contents)
    }

    pub fn to_system(&self, coordinates: &[Vector3<f64>]) -> Result<System, String> {
        if coordinates.len() != self.atoms.len() {
            return Err(format!(
                "coordinate count mismatch: got {}, expected {}",
                coordinates.len(),
                self.atoms.len()
            ));
        }

        let mut particles = Vec::with_capacity(self.atoms.len());
        for (idx, atom) in self.atoms.iter().enumerate() {
            let atom_type = self
                .atom_types
                .get(&atom.type_name)
                .ok_or_else(|| format!("missing atom type '{}'", atom.type_name))?;

            let (sigma, epsilon) = infer_lj_parameters(atom_type)?;
            let mass = atom.mass.unwrap_or(atom_type.mass);

            particles.push(Particle {
                id: atom.index,
                position: coordinates[idx],
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                lj_parameters: LJParameters {
                    epsilon,
                    sigma,
                    number_of_atoms: 1,
                },
                mass,
                energy: 0.0,
                atom_type: idx as f64,
                charge: atom.charge,
            });
        }

        let bonds = self
            .bonds
            .iter()
            .map(|b| Bond {
                atom1: b.atom1 - 1,
                atom2: b.atom2 - 1,
                k: b.force_constant,
                r0: b.length,
            })
            .collect();

        let angles = self
            .angles
            .iter()
            .map(|a| Angle {
                atom1: a.atom1 - 1,
                atom2: a.atom2 - 1,
                atom3: a.atom3 - 1,
                k: a.force_constant,
                theta0: a.theta0_deg.to_radians(),
            })
            .collect();

        let dihedrals = self
            .dihedrals
            .iter()
            .map(|d| Dihedral {
                atom1: d.atom1 - 1,
                atom2: d.atom2 - 1,
                atom3: d.atom3 - 1,
                atom4: d.atom4 - 1,
                k: d.force_constant,
                multiplicity: d.multiplicity,
                phase: d.phase_deg.to_radians(),
            })
            .collect();

        Ok(System {
            atoms: particles,
            bonds,
            angles,
            dihedrals,
            impropers: Vec::new(),
        })
    }
}

fn parse_atomtype(tokens: &[&str]) -> Result<MartiniAtomType, String> {
    if tokens.len() < 5 {
        return Err("atomtypes row requires at least 5 columns".to_string());
    }

    let name = tokens[0].to_string();
    let mass = parse_f64(tokens, 1, "atomtype mass")?;
    let charge = parse_f64(tokens, 2, "atomtype charge")?;

    // Martini files can encode either sigma/epsilon or C6/C12 in the last 2 columns.
    let maybe_a = tokens
        .get(tokens.len().saturating_sub(2))
        .ok_or_else(|| "atomtype missing nonbonded columns".to_string())?
        .parse::<f64>()
        .map_err(|e| format!("invalid atomtype nonbonded value: {e}"))?;

    let maybe_b = tokens
        .get(tokens.len().saturating_sub(1))
        .ok_or_else(|| "atomtype missing nonbonded columns".to_string())?
        .parse::<f64>()
        .map_err(|e| format!("invalid atomtype nonbonded value: {e}"))?;

    // Heuristic: sigma is usually ~0.3-0.8 nm, epsilon few kJ/mol; C12 is tiny.
    let (sigma, epsilon, c6, c12) = if maybe_a > 0.0 && maybe_a < 2.0 && maybe_b < 20.0 {
        (Some(maybe_a), Some(maybe_b), None, None)
    } else {
        (None, None, Some(maybe_a), Some(maybe_b))
    };

    Ok(MartiniAtomType {
        name,
        mass,
        charge,
        sigma,
        epsilon,
        c6,
        c12,
    })
}

fn parse_atom(tokens: &[&str]) -> Result<MartiniAtom, String> {
    if tokens.len() < 7 {
        return Err("atoms row requires at least 7 columns".to_string());
    }

    let index = parse_usize(tokens, 0, "atom index")?;
    let type_name = tokens[1].to_string();
    let charge = parse_f64(tokens, 6, "atom charge")?;
    let mass = if tokens.len() > 7 {
        Some(parse_f64(tokens, 7, "atom mass")?)
    } else {
        None
    };

    Ok(MartiniAtom {
        index,
        type_name,
        charge,
        mass,
    })
}

fn parse_bond(tokens: &[&str]) -> Result<MartiniBond, String> {
    if tokens.len() < 5 {
        return Err("bonds row requires at least 5 columns".to_string());
    }

    Ok(MartiniBond {
        atom1: parse_usize(tokens, 0, "bond atom1")?,
        atom2: parse_usize(tokens, 1, "bond atom2")?,
        length: parse_f64(tokens, 3, "bond length")?,
        force_constant: parse_f64(tokens, 4, "bond force constant")?,
    })
}

fn parse_angle(tokens: &[&str]) -> Result<MartiniAngle, String> {
    if tokens.len() < 6 {
        return Err("angles row requires at least 6 columns".to_string());
    }

    Ok(MartiniAngle {
        atom1: parse_usize(tokens, 0, "angle atom1")?,
        atom2: parse_usize(tokens, 1, "angle atom2")?,
        atom3: parse_usize(tokens, 2, "angle atom3")?,
        theta0_deg: parse_f64(tokens, 4, "angle theta0")?,
        force_constant: parse_f64(tokens, 5, "angle force constant")?,
    })
}

fn parse_dihedral(tokens: &[&str]) -> Result<MartiniDihedral, String> {
    if tokens.len() < 8 {
        return Err("dihedrals row requires at least 8 columns".to_string());
    }

    Ok(MartiniDihedral {
        atom1: parse_usize(tokens, 0, "dihedral atom1")?,
        atom2: parse_usize(tokens, 1, "dihedral atom2")?,
        atom3: parse_usize(tokens, 2, "dihedral atom3")?,
        atom4: parse_usize(tokens, 3, "dihedral atom4")?,
        phase_deg: parse_f64(tokens, 5, "dihedral phase")?,
        force_constant: parse_f64(tokens, 6, "dihedral force constant")?,
        multiplicity: parse_usize(tokens, 7, "dihedral multiplicity")?,
    })
}

fn infer_lj_parameters(atom_type: &MartiniAtomType) -> Result<(f64, f64), String> {
    if let (Some(sigma), Some(epsilon)) = (atom_type.sigma, atom_type.epsilon) {
        return Ok((sigma, epsilon));
    }

    if let (Some(c6), Some(c12)) = (atom_type.c6, atom_type.c12) {
        if c6 <= 0.0 || c12 <= 0.0 {
            return Err(format!(
                "atom type '{}' has non-positive C6/C12",
                atom_type.name
            ));
        }
        let sigma = (c12 / c6).powf(1.0 / 6.0);
        let epsilon = (c6 * c6) / (4.0 * c12);
        return Ok((sigma, epsilon));
    }

    Err(format!(
        "atom type '{}' is missing LJ parameters",
        atom_type.name
    ))
}

fn parse_f64(tokens: &[&str], index: usize, label: &str) -> Result<f64, String> {
    tokens
        .get(index)
        .ok_or_else(|| format!("missing {label}"))?
        .parse::<f64>()
        .map_err(|e| format!("failed to parse {label}: {e}"))
}

fn parse_usize(tokens: &[&str], index: usize, label: &str) -> Result<usize, String> {
    tokens
        .get(index)
        .ok_or_else(|| format!("missing {label}"))?
        .parse::<usize>()
        .map_err(|e| format!("failed to parse {label}: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_martini_itp_and_builds_system() {
        let itp = r#"
[ moleculetype ]
; name nrexcl
WAT 1

[ atomtypes ]
;name mass charge ptype sigma epsilon
P4   72.0 0.0 A 0.47 5.0

[ atoms ]
;nr type resnr residue atom cgnr charge mass
1 P4 1 WAT W1 1 0.0 72.0
2 P4 1 WAT W2 1 0.0 72.0

[ bonds ]
1 2 1 0.30 1250

[ angles ]
1 2 1 2 120.0 25.0

[ dihedrals ]
1 2 1 2 1 180.0 5.0 1
"#;

        let ff = MartiniForceField::parse_str(itp).expect("martini parsing should succeed");
        assert_eq!(ff.atoms.len(), 2);
        assert_eq!(ff.bonds.len(), 1);

        let coords = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.3, 0.0, 0.0)];
        let system = ff.to_system(&coords).expect("system build should succeed");

        assert_eq!(system.atoms.len(), 2);
        assert_eq!(system.bonds[0].atom1, 0);
        assert!((system.atoms[0].lj_parameters.sigma - 0.47).abs() < 1e-12);
        assert!((system.atoms[0].lj_parameters.epsilon - 5.0).abs() < 1e-12);
    }

    #[test]
    fn parses_c6_c12_atomtypes() {
        let itp = r#"
[ atomtypes ]
C1   72.0 0.0 A 2.0 1.0
[ atoms ]
1 C1 1 TST C1 1 0.0 72.0
"#;

        let ff = MartiniForceField::parse_str(itp).expect("martini parsing should succeed");
        let atom_type = ff.atom_types.get("C1").expect("atomtype should exist");
        let (sigma, epsilon) = infer_lj_parameters(atom_type).expect("lj inference should succeed");

        assert!(sigma > 0.0);
        assert!(epsilon > 0.0);
    }
}
