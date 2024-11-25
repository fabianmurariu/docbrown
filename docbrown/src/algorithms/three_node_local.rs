use std::collections::{HashMap};

use crate::db::{
    graph::Graph,
    program::{GlobalEvalState, LocalState, Program},
    view_api::*, edge::EdgeView,
};
use crate::algorithms::three_node_motifs::*;
use rustc_hash::FxHashMap;

pub fn star_motif_count<G:GraphViewOps>(graph:&G, v:u64, delta:i64) -> [usize;24] {
    if let Some(vertex) = graph.vertex(v) {
        let neigh_map : HashMap<u64,usize> = vertex.neighbours().iter().enumerate().map(|(num,nb)| (nb.id(), num) ).into_iter().collect();
        let exploded_edges = vertex.edges()
        .explode()
        .map(|edge| if edge.src().id()==v {star_event(neigh_map[&edge.dst().id()],1,edge.time().unwrap())} else {star_event(neigh_map[&edge.src().id()],0,edge.time().unwrap())})
        .collect::<Vec<StarEvent>>();
        let mut star_count = init_star_count(neigh_map.len());
        star_count.execute(&exploded_edges, 10);
        star_count.return_counts()
    }
    else {[0;24]}
}

pub fn twonode_motif_count<G:GraphViewOps>(graph:&G, v:u64, delta:i64) -> [u64;8] {
    let mut counts = [0;8];
    if let Some(vertex) = graph.vertex(v) {
        for nb in vertex.neighbours().iter() {
            let nb_id = nb.id();
            println!("edge {} {}",v, nb_id);
            let out = graph.edge(vertex.id(), nb_id, None);
            let inc = graph.edge(nb_id, vertex.id(),None);
            let mut all_exploded = match (out,inc) {
                (Some(o),Some(i)) => o.explode()
                .chain(i.explode())
                .map(|e| two_node_event(if e.src().id()==v {1} else {0}, e.time().unwrap()))
                .collect::<Vec<TwoNodeEvent>>(),
                (Some(o), None) => o.explode()
                .map(|e| two_node_event(1, e.time().unwrap()))
                .collect::<Vec<TwoNodeEvent>>(),
                (None, Some(i)) => i.explode()
                .map(|e| two_node_event(0, e.time().unwrap()))
                .collect::<Vec<TwoNodeEvent>>(),
                (None, None) => Vec::new()
            };
            all_exploded.sort_by_key(|e| e.time);
            for e in &all_exploded {
                println!("{} {}",e.dir, e.time);
            }
            let mut two_node_counter = init_two_node_count();
            two_node_counter.execute(&all_exploded, delta);
            let two_node_result = two_node_counter.return_counts();
            for i in 0..8 {
                counts[i]+=two_node_result[i];
            }
        }
    }
    counts
}

#[cfg(test)]
mod local_motif_test {
    use crate::db::graph::Graph;
    use crate::algorithms::three_node_local::*;

    #[test]
    fn test_init() {
        let graph = Graph::new(1);

        let vs = vec![
            (1, 2, 0),
            (2, 1, 1),
            (1, 4, 2),
            (2, 3, 3),
            (3, 2, 4),
            (3, 1, 5),
            (4, 3, 6),
            (4, 1, 7),
            (1, 5, 8),
            (3, 1, 9),
            (1, 2, 10)
        ];
        
        for (src, dst, time) in &vs {
            graph.add_edge(*time, *src, *dst, &vec![], None);
        }

        let counts = star_motif_count(&graph, 1, 100);
        print!("{:?}",counts)
    }
}