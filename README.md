# Algebraic Shifting Benchmarks

This repository contains the neccesary scripts, build instructions and examples for benchmarking the implementation of Exterior Algebraic Shifting available in Oscar.jl.

# Setup

## Setting a memory limit
Setting a ulimit before starting julia and running the script is neccessary.
Some examples were chosen due to their large memory consumption, and without setting
a limit the script will not finish.

The following command will set a ulimit of 64 gb.

```
ulimit -v 64000000
```

## Building julia environment

Using the provided Manifest.toml requires julia version 1.10.7.
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

from the project root or running the following command from a julia repl with the project environment activated.

```
    julia> include("run_times.jl")
```

## Bipartite Graphs

The examples used in the bipartite graph experiments can be found in `run_times_bipartite.jl`.

Running just the bipartite graph benchmark script can be run one after the other by either running 

```
    julia --project=. run_times_bipartite.jl
```

from the project root or running the following command from a julia repl with the project environment activated.

```
    julia> include("run_times_bipartite.jl")
```

## Surfaces

The examples used in the surfaces experiments can be found in the `examples/surfaces` folder.
These simplicial complexes were taken from Frank Lutz's webpage, for a complete list of available triangualtions of surfaces and more see [here](https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/)

Running just the surfaces benchmark script can be run by either the following command from the project root 

```
    julia --project=. run_times_surfaces.jl
```

or running 

```
    julia> include("run_times_surfaces.jl")
```

from a julia repl with the project environment activated.

## Non surfaces

The examples used in the surfaces experiments can be found in the `examples/non_surfaces` folder.

Running just the non surfaces benchmark script can be run by either the following command from the project root 

```
    julia --project=. run_times_non_surfaces.jl
```

or running 

```
    julia> include("run_times_non_surfaces.jl")
```

from a julia repl with the project environment activated.
