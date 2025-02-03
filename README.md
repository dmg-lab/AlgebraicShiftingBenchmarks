# Algebraic Shifting Benchmarks

This repository contains the neccesary scripts and examples for benchmarking the implementation of Exterior Algebraic Shifting available in Oscar.jl.

This is the software companion to the article "Faster Algebraic Shifting" by Antony Della Vecchia, Michael Joswig and Fabian Lenzen.

We also include the file `examples.m2`, which was used to run a comparison with the Monte Carlo algorithm implemented in https://github.com/ank1494/ext-shifting .

# Setup

## Setting a memory limit
Setting a ulimit before starting julia and running the script is neccessary.
Some examples were chosen due to their large memory consumption, and without setting
a limit the script will not finish.

The following command will set a ulimit of 80GiB and disables swapping:

```
    ulimit -v 80000000 -x 0
```

## Building julia environment

Using the provided Manifest.toml requires julia version 1.10.8.
To start julia from the root of the project and activate the environment for the project run.

```
    julia --project=.
```

Running the following command from the julia repl will setup the necessary dependancies. 
We have included a Manifest.toml so that a particular branch of Oscar is being used for benchmarking. 

```
    julia> ]up
```

# Experiments

Each experiment has it's own script, where each script is more or less structured the same.
Running each script will output a csv file, named according to the experiment and machine it was run on.
All scripts can be run one after the other by either running 

```
    julia --project=. run_times.jl
```

We include the following experiments:

 * `run_times_bipartite.jl`: different bipartite graphs
 * `run_times_surfaces.jl`: surface triangulations (in the `examples/surfaces`).
   These were taken from Frank Lutz's [compilation of manifold triangulations](https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/).
 * `run_times_surfaces_extra.jl`: A few surface triangulations are made subject to more experiments with more varying parameters.
 * `run_times_non_surfaces.hl` manually built simplicial complexes with prescribed H1 (in `examples/non_surfaces`).
   See article for details how these are built.