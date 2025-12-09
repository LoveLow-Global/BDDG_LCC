use anyhow::Result;
use rand::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use serde::Serialize;
use std::cmp::max;
use std::fs::OpenOptions;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

// Raw Graph Structure
// We store both directions to avoid re-building them later.
struct RawGraph {
    n: usize,
    fwd: Vec<Vec<u32>>, // Forward edges
    bwd: Vec<Vec<u32>>, // Backward edges
}

impl RawGraph {
    fn new(n: usize) -> Self {
        // Pre-allocate the outer vectors
        Self {
            n,
            fwd: vec![Vec::new(); n],
            bwd: vec![Vec::new(); n],
        }
    }

    // Extract the subgraph induced by a set of nodes
    // Returns a new RawGraph with remapped indices 0..subset_size
    fn induced_subgraph(&self, nodes: &[u32]) -> RawGraph {
        let sub_n = nodes.len();
        let mut sub_g = RawGraph::new(sub_n);

        // Map original_id -> new_id
        // If graph is not dense, HashMap will be better. To be edited.
        let mut node_map = vec![u32::MAX; self.n];
        for (new_idx, &orig_node) in nodes.iter().enumerate() {
            node_map[orig_node as usize] = new_idx as u32;
        }

        for (new_u, &orig_u) in nodes.iter().enumerate() {
            let u_idx = new_u as u32;
            // Iterate original edges
            for &orig_v in &self.fwd[orig_u as usize] {
                // If neighbor is in the subgraph, add edge
                let v_idx = node_map[orig_v as usize];
                if v_idx != u32::MAX {
                    sub_g.fwd[new_u].push(v_idx);
                    sub_g.bwd[v_idx as usize].push(u_idx);
                }
            }
        }
        sub_g
    }
}

#[derive(Serialize)]
struct Record {
    n: usize,
    p: f64,
    rep: usize,
    diameter: u32,
}

/// Generate graph directly into Adjacency Lists without Petgraphs
fn ber_directed_divisor_graph(n: usize, p: f64, rng: &mut StdRng) -> RawGraph {
    let mut g = RawGraph::new(n);

    // Capacity reservation to reduce re-allocations
    // Avg degree is approx 2 * ln(n)
    let est_cap = 2*(p * (n as f64).ln()) as usize;
    for i in 0..n {
        g.fwd[i].reserve(est_cap);
        g.bwd[i].reserve(est_cap);
    }

    for i_val in 1..=n {
        let u = (i_val - 1) as u32;
        
        let mut j_val = 2 * i_val;
        while j_val <= n {
            let v = (j_val - 1) as u32;

            if rng.random_bool(p) {
                // u -> v
                g.fwd[u as usize].push(v);
                g.bwd[v as usize].push(u);
            } else {
                // v -> u
                g.fwd[v as usize].push(u);
                g.bwd[u as usize].push(v);
            }
            j_val += i_val;
        }
    }
    g
}

// Iterative Kosaraju
fn find_largest_scc(g: &RawGraph) -> Vec<u32> {
    let n = g.n;
    let mut visited = vec![false; n];
    let mut stack = Vec::with_capacity(n); // Finish order
    
    // DFS on Forward Graph to compute finish order
    // Iterative DFS
    let mut dfs_stack = Vec::with_capacity(n); // (node, child_idx)
    
    for start_node in 0..n {
        if visited[start_node] { continue; }
        
        // Push start node
        dfs_stack.push((start_node as u32, 0));
        visited[start_node] = true;

        while let Some((u, child_idx)) = dfs_stack.last_mut() {
            let u_usize = *u as usize;
            let edges = &g.fwd[u_usize];
            
            if *child_idx < edges.len() {
                let v = edges[*child_idx];
                *child_idx += 1;
                let v_usize = v as usize;
                if !visited[v_usize] {
                    visited[v_usize] = true;
                    dfs_stack.push((v, 0));
                }
            } else {
                // All children visited, finish u
                stack.push(*u);
                dfs_stack.pop();
            }
        }
    }

    //  DFS on Backward Graph using finish order
    visited.fill(false);
    let mut largest_component = Vec::new();
    let mut current_component = Vec::new();

    while let Some(root) = stack.pop() {
        if visited[root as usize] { continue; }
        
        current_component.clear();
        // Start DFS from root on BWD graph
        dfs_stack.push((root, 0));
        visited[root as usize] = true;
        current_component.push(root);

        while let Some((u, child_idx)) = dfs_stack.last_mut() {
            let u_usize = *u as usize;
            let edges = &g.bwd[u_usize]; // Using Backward edges
            
            if *child_idx < edges.len() {
                let v = edges[*child_idx];
                *child_idx += 1;
                let v_usize = v as usize;
                if !visited[v_usize] {
                    visited[v_usize] = true;
                    current_component.push(v);
                    dfs_stack.push((v, 0));
                }
            } else {
                dfs_stack.pop();
            }
        }

        if current_component.len() > largest_component.len() {
            largest_component.clear();
            largest_component.extend_from_slice(&current_component);
        }
    }

    largest_component
}


// Diameter BFS
fn diameter_exact(g: &RawGraph) -> u32 {
    let n = g.n;
    if n == 0 { return 0; }
    
    let mut active = vec![true; n];
    let mut visited = vec![false; n];
    let mut queue = vec![0u32; n];
    let mut distbuf = vec![-1i32; n];
    let mut diam = 0;

    let mut vs: Vec<usize> = (0..n).collect();
    // Sort by degree
    vs.sort_unstable_by_key(|&v| {
        let deg = g.fwd[v].len() + g.bwd[v].len();
        std::cmp::Reverse(deg)
    });

    for &u in &vs {
        if !active[u] { continue; }

        // Forward BFS
        visited.fill(false);
        unsafe { *visited.get_unchecked_mut(u) = true; }
        queue[0] = u as u32;
        
        let mut front = 0;
        let mut back = 1;
        let mut level_end = 1;
        let mut e = 0;

        while front < back {

            let v = unsafe { *queue.get_unchecked(front) } as usize;
            front += 1;

            let neighbors = unsafe { g.fwd.get_unchecked(v) };
            for &w in neighbors {
                let w_idx = w as usize;
                // We know w_idx < n from construction
                let vis = unsafe { visited.get_unchecked_mut(w_idx) };
                if !*vis {
                    *vis = true;
                    unsafe { *queue.get_unchecked_mut(back) = w };
                    back += 1;
                }
            }

            if front > level_end && front < back {
                e += 1;
                level_end = back - 1;
            }
        }
        
        diam = max(diam, e);

        // Backward BFS : Pruning
        let dmax = diam - e;
        if dmax >= 0 {
            distbuf.fill(-1);
            unsafe { *distbuf.get_unchecked_mut(u) = 0; }
            queue[0] = u as u32;
            front = 0;
            back = 1;

            while front < back {
                let v = unsafe { *queue.get_unchecked(front) } as usize;
                front += 1;
                
                let d_v = unsafe { *distbuf.get_unchecked(v) };
                if d_v >= dmax { continue; }

                let neighbors = unsafe { g.bwd.get_unchecked(v) };
                for &w in neighbors {
                    let w_idx = w as usize;
                    let d_w = unsafe { distbuf.get_unchecked_mut(w_idx) };
                    if *d_w < 0 {
                        *d_w = d_v + 1;
                        unsafe { *queue.get_unchecked_mut(back) = w };
                        back += 1;
                    }
                }
            }

            for v in 0..n {
                // If active[v] is true
                if unsafe { *active.get_unchecked(v) } {
                     let d = unsafe { *distbuf.get_unchecked(v) };
                     if d >= 0 && (d + e) <= diam {
                         unsafe { *active.get_unchecked_mut(v) = false; }
                     }
                }
            }
        }

        if !active.iter().any(|&x| x) { break; }
    }

    diam as u32
}

fn main() -> Result<()> {
    let n = 1 << 10; 
    let ps: Vec<f64> = (1..=50).map(|x| x as f64 / 100.0).collect();
    let n_replicate = 10;

    let filename = format!("4_2to{}_diameter_results_rs.csv", (n as f64).log2() as u32);
    let path = Path::new(&filename);
    let file_exists = path.exists();
    let file = OpenOptions::new().append(true).create(true).open(path)?;
    
    if !file_exists {
        let mut wtr = csv::Writer::from_writer(&file);
        wtr.write_record(&["n", "p", "rep", "diameter"])?;
        wtr.flush()?;
    }

    let file_mutex = Arc::new(Mutex::new(file));

    println!("Starting simulation with N={}, threads={}", n, rayon::current_num_threads());
    let start_total = Instant::now();

    for p in ps {
        let loop_start = Instant::now();
        
        let results: Vec<Record> = (1..=n_replicate).into_par_iter().map(|rep| {
            let mut rng = StdRng::seed_from_u64(42 + (p * 100.0) as u64 + rep as u64);

            // Raw Graph
            let g = ber_directed_divisor_graph(n, p, &mut rng);

            // Iterative SCC
            let largest_nodes = find_largest_scc(&g);

            if !largest_nodes.is_empty() {
                // Induce Subgraph
                let subg = g.induced_subgraph(&largest_nodes);
                
                // Exact Diameter
                let d = diameter_exact(&subg);
                Record { n, p, rep, diameter: d }
            } else {
                Record { n, p, rep, diameter: 0 }
            }

        }).collect();

        {
            let file_guard = file_mutex.lock().unwrap();
            let mut wtr = csv::WriterBuilder::new().has_headers(false).from_writer(&*file_guard);
            for r in &results {
                wtr.serialize(r).expect("CSV serialize error");
            }
            wtr.flush().expect("CSV flush error");
        }

        let duration = loop_start.elapsed();
        println!("Completed p = {:.2} with {} reps in {:.2?}", p, n_replicate, duration);
    }

    println!("Total time: {:.4?}", start_total.elapsed());
    Ok(())
}