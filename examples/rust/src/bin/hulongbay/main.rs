#![allow(unused_imports)]
use std::collections::HashMap;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;
use std::{env, thread};

use chrono::{DateTime, Utc};
use docbrown_core::tgraph::TemporalGraph;
use docbrown_core::utils;
use docbrown_core::{Direction, Prop};
use docbrown_db::csv_loader::csv::CsvLoader;
use regex::Regex;
use serde::Deserialize;
use std::fs::File;
use std::io::{prelude::*, BufReader, LineWriter};
use std::time::Instant;

use docbrown_db::graph::Graph;
use docbrown_db::view_api::*;

#[derive(Deserialize, std::fmt::Debug)]
pub struct Edge {
    _unknown0: i64,
    _unknown1: i64,
    _unknown2: i64,
    src: u64,
    dst: u64,
    time: i64,
    _unknown3: u64,
    amount_usd: u64,
}

#[derive(Debug)]
pub struct MissingArgumentError;

impl Display for MissingArgumentError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Failed to provide the path to the hulongbay data directory"
        )
    }
}

impl Error for MissingArgumentError {}

#[derive(Debug)]
pub struct GraphEmptyError;

impl Display for GraphEmptyError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "The graph was empty but data was expected.")
    }
}

impl Error for GraphEmptyError {}

pub fn loader(data_dir: &Path) -> Result<Graph, Box<dyn Error>> {
    let encoded_data_dir = data_dir.join("graphdb.bincode");
    if encoded_data_dir.exists() {
        let now = Instant::now();
        let g = Graph::load_from_file(encoded_data_dir.as_path())?;

        println!(
            "Loaded graph from path {} with {} vertices, {} edges, took {} seconds",
            encoded_data_dir.display(),
            g.num_vertices(),
            g.num_edges(),
            now.elapsed().as_secs()
        );

        Ok(g)
    } else {
        let g = Graph::new(16);

        let now = Instant::now();

        CsvLoader::new(data_dir).load_into_graph(&g, |sent: Edge, g: &Graph| {
            let src = sent.src;
            let dst = sent.dst;
            let time = sent.time;

            g.add_edge(
                time,
                src,
                dst,
                &vec![("amount".to_owned(), Prop::U64(sent.amount_usd))],
            )
        })?;

        println!(
            "Loaded graph from CSV data files {} with {} vertices, {} edges which took {} seconds",
            encoded_data_dir.display(),
            g.num_vertices(),
            g.num_edges(),
            now.elapsed().as_secs()
        );

        g.save_to_file(encoded_data_dir)?;
        Ok(g)
    }
}

fn try_main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let data_dir = Path::new(args.get(1).ok_or(MissingArgumentError)?);

    let graph = loader(data_dir)?;

    let now = Instant::now();
    let num_edges: usize = graph.vertices().map(|v| v.out_degree()).sum();
    println!(
        "Counting edges by summing degrees returned {} in {} seconds",
        num_edges,
        now.elapsed().as_secs()
    );
    let earliest_time = graph.earliest_time().ok_or(GraphEmptyError)?;
    let latest_time = graph.latest_time().ok_or(GraphEmptyError)?;
    println!("graph time range: {}-{}", earliest_time, latest_time);
    Ok(())
}

fn main() {
    if let Err(e) = try_main() {
        eprintln!("Failed: {}", e);
        std::process::exit(1)
    }
}
