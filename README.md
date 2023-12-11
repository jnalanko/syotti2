# Syotti 2 Early Access Beta

# Installing

First, install the [Rust](https://www.rust-lang.org/tools/install). Then:

```
git submodule update --init
cargo build --release
```

This produces the binary to `target/release/unitig_flipper`. If you want to install the program to `$PATH`, run `cargo install --path .`

# Usage

```
Usage: syotti2 [OPTIONS] <COMMAND>

Commands:
  design    
  coverage  
  help      Print this message or the help of the given subcommand(s)

Options:
  -j, --threads <threads>  Number of threads to use [default: 8]
  -h, --help               Print hel
```

## Subcommand `design`

```
Usage: syotti2 design [OPTIONS] --targets <targets> --output <output>

Options:
  -j, --threads <threads>
          Number of threads to use [default: 8]
  -t, --targets <targets>
          Input FASTA or FASTQ file, possibly gzipped
  -o, --output <output>
          FASTA outputfile for the baits
  -L, --bait-len <bait-length>
          Bait length [default: 120]
  -d, --hamming-distance <hamming-distance>
          Maximum number of mismatches allowed in matching [default: 40]
  -g, --seed-len <seed-len>
          Length of seeds in matching [default: 20]
  -m, --minimizer-len <minimizer-len>
          Length of minimizers in indexing. Must be less or equal to seed-len [default: 12]
  -c, --cutoff <cutoff>
          Stop the algorithm when this coverage fraction is reached [default: 1.0]
  -r, --randomize
          Randomize the processing order in the greedy algorithm
  -h, --help
          Print help
```

## Subcommand `coverage`

```
Usage: syotti2 coverage [OPTIONS] --baits <baits> --targets <targets> --output <output>

Options:
  -b, --baits <baits>
          FASTA file containing the baits
  -j, --threads <threads>
          Number of threads to use [default: 8]
  -t, --targets <targets>
          FASTA or FASTQ file, possibly gzipped
  -r, --resolution <resolution>
          Computes a moving average over the coverage vectors and samples those at evenly spaced intervals to get this many data points per sequence.
  -d, --hamming-distance <hamming-distance>
          Maximum number of mismatches allowed in matching [default: 40]
  -g, --seed-len <seed-len>
          Length of seeds in matching [default: 20]
  -m, --minimizer-len <minimizer-len>
          Length of minimizers in indexing. Must be less or equal to seed-len [default: 12]
  -o, --output <output>
          Output csv file for the coverage data.
      --mismatch-out <mismatch-out>
          Output csv file for the mismatch data.
  -h, --help
          Print help (see more with '--help')
```


