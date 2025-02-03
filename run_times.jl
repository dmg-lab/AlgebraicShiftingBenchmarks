# Run this file with (for example):
#   ulimit -v 80000000 -x 0
#   julia --heap-size-hint=1 --project=. run_times.jl

include("run_times_bipartite.jl")
include("run_times_surfaces.jl")
include("run_times_surfaces_extra.jl")
include("run_times_non_surfaces.jl")

