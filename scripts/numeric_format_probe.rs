//! Standalone probe of the exact Rust Display formatting used by plugin_cache/csq.rs.
//! Compile with rustc -O; stdin is value<TAB>frequency, argv[1] is fixed scale.
use std::io::{self, BufRead};

fn main() {
    let scale: usize = std::env::args().nth(1).unwrap().parse().unwrap();
    let mut counts = [0u64; 6];
    let mut examples: [Vec<String>; 5] = Default::default();
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let (s, count) = line.rsplit_once('\t').unwrap();
        let count: u64 = count.parse().unwrap();
        counts[0] += count;
        let Ok(v64) = s.parse::<f64>() else {
            counts[1] += count;
            if examples[0].len() < 5 {
                examples[0].push(s.to_owned());
            }
            continue;
        };
        let v32 = s.parse::<f32>().unwrap();
        let outputs = [
            v32.to_string(),
            v64.to_string(),
            fixed32(v32, scale),
            fixed64(v64, scale),
        ];
        for (i, rendered) in outputs.into_iter().enumerate() {
            if rendered != s {
                counts[i + 2] += count;
                if examples[i + 1].len() < 5 {
                    examples[i + 1].push(format!("{s} -> {rendered}"));
                }
            }
        }
    }
    println!("{}", counts.map(|x| x.to_string()).join("\t"));
    for (i, items) in examples.iter().enumerate() {
        for item in items {
            println!("{i}\t{item}");
        }
    }
}

// Experimental renderer policies: 100 = CADD PHRED; 101 = GERP source spelling.
fn phred_scale(v: f64) -> usize {
    if v < 10.0 {
        3
    } else if v < 20.0 {
        2
    } else if v < 30.0 {
        1
    } else {
        0
    }
}
fn gerp(decimal: String, scientific: String, small: bool) -> String {
    if small {
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let mantissa = if mantissa.contains('.') {
            mantissa.to_owned()
        } else {
            format!("{mantissa}.0")
        };
        format!("{mantissa}E{exponent}")
    } else if decimal.contains('.') {
        decimal
    } else {
        format!("{decimal}.0")
    }
}
fn fixed32(v: f32, scale: usize) -> String {
    if scale == 101 {
        return gerp(v.to_string(), format!("{v:e}"), v != 0.0 && v.abs() < 0.001);
    }
    let scale = if scale == 100 {
        phred_scale(v as f64)
    } else {
        scale
    };
    format!("{v:.scale$}")
}
fn fixed64(v: f64, scale: usize) -> String {
    if scale == 101 {
        return gerp(v.to_string(), format!("{v:e}"), v != 0.0 && v.abs() < 0.001);
    }
    let scale = if scale == 100 { phred_scale(v) } else { scale };
    format!("{v:.scale$}")
}
