use serde::{Deserialize, Serialize};

#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;
extern crate core;

mod adj;
pub mod agg;
mod bitset;
mod edge_layer;
mod lazy_vec;
pub mod lsm;
mod misc;
mod props;
mod sorted_vec_map;
pub mod state;
pub mod tadjset;
mod tcell;
pub mod tgraph;
pub mod tgraph_shard;
mod tprop;
pub mod utils;
pub mod vertex;

// Denotes edge direction
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Direction {
    OUT,
    IN,
    BOTH,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Clone)]
pub enum Prop {
    Str(String),
    I32(i32),
    I64(i64),
    U32(u32),
    U64(u64),
    F32(f32),
    F64(f64),
    Bool(bool),
}
