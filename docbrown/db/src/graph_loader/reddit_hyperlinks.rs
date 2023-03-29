//! Load the Reddit hyperlinks dataset into a graph.
//! The dataset is available at https://snap.stanford.edu/data/soc-redditHyperlinks-body.tsv
//! The hyperlink network represents the directed connections between two subreddits (a subreddit
//! is a community on Reddit). We also provide subreddit embeddings. The network is extracted
//! from publicly available Reddit data of 2.5 years from Jan 2014 to April 2017.
//! *NOTE: It may take a while to download the dataset
//!
//! ## Dataset statistics
//! * Number of nodes (subreddits) 	55,863
//! * Number of edges (hyperlink between subreddits) 	858,490
//! * Edge weights (label of hyperlink) 	-1 or +1
//! * Edge attributes 	Text property vectors
//! * Timespan 	Jan 2014 - April 2017
//!
//! ## Source
//! S. Kumar, W.L. Hamilton, J. Leskovec, D. Jurafsky. Community Interaction and Conflict
//! on the Web. World Wide Web Conference, 2018.
//!
//! ## Properties
//!
//!  * SOURCE_SUBREDDIT: the subreddit where the link originates
//!  * TARGET_SUBREDDIT: the subreddit where the link ends
//!  * POST_ID: the post in the source subreddit that starts the link
//!  * TIMESTAMP: time time of the post
//!  * POST_LABEL: label indicating if the source post is explicitly negative towards the target
//! post. The value is -1 if the source is negative towards the target, and 1 if it is neutral or
//! positive. The label is created using crowd-sourcing and training a text based classifier, and
//! is better than simple sentiment analysis of the posts. Please see the reference paper for details.
//!  * POST_PROPERTIES: a vector representing the text properties of the source post, listed as a
//! list of comma separated numbers. This can be found on the source website
//!
//! Example:
//! ```rust
//! use docbrown_db::graph_loader::reddit_hyperlinks::reddit_graph;
//! use docbrown_db::graph::Graph;
//! use docbrown_db::view_api::*;
//!
//! let graph = reddit_graph(1, 120);
//!
//! println!("The graph has {:?} vertices", graph.num_vertices());
//! println!("The graph has {:?} edges", graph.num_edges());
//! ```

use crate::{graph::Graph, graph_loader::fetch_file};
use chrono::*;
use docbrown_core::Prop;
use itertools::Itertools;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::path::PathBuf;

/// Download the dataset and return the path to the file
pub fn reddit_file(timeout: u64) -> Result<PathBuf, Box<dyn std::error::Error>> {
    fetch_file(
        "reddit.tsv",
        "https://snap.stanford.edu/data/soc-redditHyperlinks-body.tsv",
        timeout,
    )
}

/// Read the file line by line
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// Load the Reddit hyperlinks dataset into a graph and return it
pub fn reddit_graph(shards: usize, timeout: u64) -> Graph {
    let graph = {
        let g = Graph::new(shards);

        if let Ok(path) = reddit_file(timeout) {
            if let Ok(lines) = read_lines(path.as_path()) {
                // Consumes the iterator, returns an (Optional) String
                for line in lines.dropping(1) {
                    if let Ok(reddit) = line {
                        let reddit: Vec<&str> = reddit.split("	").collect();
                        let src_id = &reddit[0];
                        let dst_id = &reddit[1];
                        let post_id = reddit[2].to_string();

                        match NaiveDateTime::parse_from_str(reddit[3], "%Y-%m-%d %H:%M:%S") {
                            Ok(time) => {
                                let time = time.timestamp();
                                let post_label: i32 = reddit[4].parse::<i32>().unwrap();
                                let post_properties: Vec<f64> = reddit[5]
                                    .split(",")
                                    .map(|s| s.parse::<f64>().unwrap())
                                    .collect();
                                let edge_properties = &vec![
                                    ("post_label".to_string(), Prop::I32(post_label)),
                                    ("post_id".to_string(), Prop::Str(post_id)),
                                    ("word_count".to_string(), Prop::F64(post_properties[7])),
                                    ("long_words".to_string(), Prop::F64(post_properties[9])),
                                    ("sentences".to_string(), Prop::F64(post_properties[13])),
                                    ("readability".to_string(), Prop::F64(post_properties[17])),
                                    (
                                        "positive_sentiment".to_string(),
                                        Prop::F64(post_properties[17]),
                                    ),
                                    (
                                        "negative_sentiment".to_string(),
                                        Prop::F64(post_properties[17]),
                                    ),
                                    (
                                        "compound_sentiment".to_string(),
                                        Prop::F64(post_properties[17]),
                                    ),
                                ];
                                g.add_vertex(time, src_id.clone(), &vec![])
                                    .map_err(|err| println!("{:?}", err))
                                    .ok();
                                g.add_vertex(time, dst_id.clone(), &vec![])
                                    .map_err(|err| println!("{:?}", err))
                                    .ok();
                                g.add_edge(time, src_id.clone(), dst_id.clone(), edge_properties);
                            }
                            Err(e) => {
                                println!("{}", e)
                            }
                        }
                    }
                }
            }
        };

        g
    };
    graph
}
// #[cfg(test)]
// mod reddit_test {
//     use crate::graph_loader::reddit_hyperlinks::reddit_graph;
//     use crate::view_api::GraphViewOps;
//
//     #[test]
//     fn check_data() {
//         if let Ok(g) = reddit_graph(1, 600) {
//             println!("{} {}",g.num_vertices(),g.num_edges());
//         };
//     }
// }
//
