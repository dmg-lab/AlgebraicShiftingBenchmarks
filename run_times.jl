#!/usr/bin/env prlimit -v 32000000000 -x 0 -- julia --heap-size-hint=1 --project=.

include("run_times_bipartite.jl")
include("run_times_surfaces.jl")
include("run_times_non_surfaces.jl")



