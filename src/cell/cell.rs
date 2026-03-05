use std::f64;

/// 3D vector helper (minimal).
#[derive(Clone, Copy, Debug)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

/// Periodic wrap into [0, L).
fn wrap_0_l(x: f64, l: f64) -> f64 {
    // Works for negative too
    let mut y = x % l;
    if y < 0.0 {
        y += l;
    }
    y
}

/// Minimum-image displacement for periodic boundary conditions.
/// Returns dx in (-L/2, L/2].
fn min_image(mut dx: f64, l: f64) -> f64 {
    let half = 0.5 * l;
    if dx > half {
        dx -= l;
    }
    if dx <= -half {
        dx += l;
    }
    dx
}

/// Cell-list structure:
/// - head[cell] = index of first particle in that cell, or None
/// - next[i] = next particle index in the same cell, or None
///
/// This is the "linked list in arrays" approach.
struct CellList {
    // domain
    box_len: Vec3,
    // cutoff and cell geometry
    cutoff: f64,
    cell_size: f64, // typically >= cutoff
    nx: usize,
    ny: usize,
    nz: usize,

    // linked-list storage
    head: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
}

impl CellList {
    /// Create a cell list. We set cell_size = cutoff (common choice),
    /// and number of cells along each axis = floor(L / cell_size).
    /// Ensure nx,ny,nz >= 1.
    fn new(box_len: Vec3, cutoff: f64) -> Self {
        let cell_size = cutoff;
        let nx = (box_len.x / cell_size).floor().max(1.0) as usize;
        let ny = (box_len.y / cell_size).floor().max(1.0) as usize;
        let nz = (box_len.z / cell_size).floor().max(1.0) as usize;

        let ncell = nx * ny * nz;

        Self {
            box_len,
            cutoff,
            cell_size,
            nx,
            ny,
            nz,
            head: vec![None; ncell],
            next: Vec::new(), // sized on rebuild
        }
    }

    #[inline]
    fn ncell(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Convert (cx,cy,cz) -> linear cell id.
    #[inline]
    fn cell_id(&self, cx: usize, cy: usize, cz: usize) -> usize {
        // row-major: x fastest
        (cz * self.ny + cy) * self.nx + cx
    }

    /// Wrap integer cell coordinate with periodic boundaries.
    #[inline]
    fn wrap_c(&self, c: isize, n: usize) -> usize {
        // periodic wrap in [0, n)
        let n = n as isize;
        let mut r = c % n;
        if r < 0 {
            r += n;
        }
        r as usize
    }

    /// Map particle position to a cell coordinate (cx,cy,cz).
    /// Assumes positions can be outside box; we wrap them.
    #[inline]
    fn pos_to_cell(&self, p: Vec3) -> (usize, usize, usize) {
        let x = wrap_0_l(p.x, self.box_len.x);
        let y = wrap_0_l(p.y, self.box_len.y);
        let z = wrap_0_l(p.z, self.box_len.z);

        let cx = (x / self.cell_size).floor() as usize % self.nx;
        let cy = (y / self.cell_size).floor() as usize % self.ny;
        let cz = (z / self.cell_size).floor() as usize % self.nz;
        (cx, cy, cz)
    }

    /// Rebuild linked list from scratch in O(N).
    ///
    /// This is the "what is wrong with rebuilding the cell list?" part:
    /// it's cheap and common to rebuild every step (or every few steps with a skin).
    fn rebuild(&mut self, positions: &[Vec3]) {
        // clear heads
        self.head.fill(None);

        // reset next pointers for each particle
        self.next.clear();
        self.next.resize(positions.len(), None);

        for (i, &p) in positions.iter().enumerate() {
            let (cx, cy, cz) = self.pos_to_cell(p);
            let c = self.cell_id(cx, cy, cz);

            // Insert particle i at the head of the cell's linked list:
            // next[i] becomes old head, head becomes i.
            self.next[i] = self.head[c];
            self.head[c] = Some(i);
        }
    }

    /// Iterate particle indices in a cell by walking the linked list.
    #[inline]
    fn iter_cell<'a>(&'a self, cell: usize) -> CellIter<'a> {
        CellIter {
            cl: self,
            cur: self.head[cell],
        }
    }

    /// Neighbor traversal: for each particle i, visit candidates j in the 27 neighboring cells.
    ///
    /// This yields pairs (i,j) with j>i (no double-counting),
    /// and you can do your distance check + force calc inside the callback.
    fn for_each_neighbor_pair<F>(&self, positions: &[Vec3], mut f: F)
    where
        F: FnMut(usize, usize, Vec3, f64), // (i, j, dr, r2)
    {
        let rc2 = self.cutoff * self.cutoff;

        // Loop over all cells
        for cz in 0..self.nz {
            for cy in 0..self.ny {
                for cx in 0..self.nx {
                    let c0 = self.cell_id(cx, cy, cz);

                    // For each particle i in this cell
                    let mut it_i = self.iter_cell(c0);
                    while let Some(i) = it_i.next() {
                        let pi = positions[i];

                        // Check 27 neighboring cells (including itself)
                        for dz in -1isize..=1 {
                            for dy in -1isize..=1 {
                                for dx in -1isize..=1 {
                                    let nx = self.wrap_c(cx as isize + dx, self.nx);
                                    let ny = self.wrap_c(cy as isize + dy, self.ny);
                                    let nz = self.wrap_c(cz as isize + dz, self.nz);
                                    let c1 = self.cell_id(nx, ny, nz);

                                    // Walk candidates j in neighbor cell
                                    let mut it_j = self.iter_cell(c1);
                                    while let Some(j) = it_j.next() {
                                        // avoid double-count and self-pair
                                        if j <= i {
                                            continue;
                                        }

                                        let pj = positions[j];

                                        // Minimum-image displacement
                                        let dr = Vec3::new(
                                            min_image(pj.x - pi.x, self.box_len.x),
                                            min_image(pj.y - pi.y, self.box_len.y),
                                            min_image(pj.z - pi.z, self.box_len.z),
                                        );
                                        let r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                                        if r2 <= rc2 {
                                            f(i, j, dr, r2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Iterator over a cell's particle indices by following `next`.
struct CellIter<'a> {
    cl: &'a CellList,
    cur: Option<usize>,
}

impl<'a> Iterator for CellIter<'a> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let i = self.cur?;
        self.cur = self.cl.next[i];
        Some(i)
    }
}

//// ----------------------------
//// Example usage
//// ----------------------------
//fn main() {
//    // Example box and cutoff (units arbitrary)
//    let box_len = Vec3::new(10.0, 10.0, 10.0);
//    let cutoff = 2.5;
//
//    let mut positions = vec![
//        Vec3::new(1.0, 1.0, 1.0),
//        Vec3::new(2.0, 1.2, 1.1),
//        Vec3::new(8.9, 9.1, 9.0), // near boundary
//        Vec3::new(0.2, 0.1, 0.3), // periodic neighbor with particle 2 potentially
//    ];
//
//    let mut cl = CellList::new(box_len, cutoff);
//
//    // Rebuild every step (or every k steps)
//    cl.rebuild(&positions);
//
//    // Count neighbor pairs (within cutoff)
//    let mut count = 0usize;
//    cl.for_each_neighbor_pair(&positions, |i, j, dr, r2| {
//        count += 1;
//        // Insert your force calc here (LJ etc.)
//        println!(
//            "pair ({i},{j}) r2={r2:.4} dr=({:.3},{:.3},{:.3})",
//            dr.x, dr.y, dr.z
//        );
//    });
//
//    println!("pairs within cutoff = {count}");
//
//    // Example "move particles" then rebuild again
//    positions[0].x += 0.05;
//    positions[1].x += 0.05;
//    cl.rebuild(&positions);
//}
//
