#!/usr/bin/env prlimit -v 32000000000 -x 0 -- julia --heap-size-hint=1 --project=.
using Revise, Oscar, DataStructures

include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
computer_name = gethostname()
bipartite_table_path = joinpath(root_of_project, "bipartite_table_test_$computer_name.csv")
surfaces_table_path = joinpath(root_of_project, "surfaces_table_test_$computer_name.csv")
non_surfaces_table_path = joinpath(root_of_project, "non_surfaces_table_test_$computer_name.csv")

# first table bipartite one field all algorithms
# las vegas needs number of computation
examples = [
  ("\$K_{32}\$", complete_bipartite_graph(3, 2)),
  ("\$K_{23}\$", complete_bipartite_graph(2, 3)),
  ("\$K_{33}\$", complete_bipartite_graph(3, 3)),
  ("\$K_{42}\$", complete_bipartite_graph(4, 2)),
  ("\$K_{24}\$", complete_bipartite_graph(2, 4)),
  ("\$K_{43}\$", complete_bipartite_graph(4, 3)),
  ("\$K_{34}\$", complete_bipartite_graph(3, 4)),
  ("\$K_{44}\$", complete_bipartite_graph(4, 4)),
  ("\$K_{52}\$", complete_bipartite_graph(5, 2)),
  ("\$K_{25}\$", complete_bipartite_graph(2, 5)),
  ("\$K_{53}\$", complete_bipartite_graph(5, 3)),
  ("\$K_{35}\$", complete_bipartite_graph(3, 5)),
  ("\$K_{54}\$", complete_bipartite_graph(5, 4)),
  ("\$K_{45}\$", complete_bipartite_graph(4, 5)),
  ("\$K_{55}\$", complete_bipartite_graph(5, 5)),
  ("\$K_{62}\$", complete_bipartite_graph(6, 2)),
  ("\$K_{26}\$", complete_bipartite_graph(2, 6)),
  ("\$K_{63}\$", complete_bipartite_graph(6, 3)),
  ("\$K_{36}\$", complete_bipartite_graph(3, 6)),
  ("\$K_{64}\$", complete_bipartite_graph(6, 4)),
  ("\$K_{46}\$", complete_bipartite_graph(4, 6)),
  ("\$K_{65}\$", complete_bipartite_graph(6, 5)),
]

# Algorithms to run, and the labels for the columns for that algorithm.
algorithms = [
  "av"   => ["avTime", "avMemory"],
  "hv"   => ["hvTime", "hvMemory"],
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  "avf"  => ["avfTime", "avfMemory"],
  "hvf"  => ["hvfTime", "hvfMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]

# If `useremote` is set, worker processes are prepared for remote execution asynchronously while the benchmark is running.
# To be safe from memory leaks etc., the workers are removed at the end of each benchmark.
# Start a few workers to have them ready for the benchmark.
useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

add_ref_labels(alg, alg_labels) = [alg_labels; alg .* ref_labels]
function bipartite_table(examples)
  open(bipartite_table_path, "w") do f
    # Print the table headers
    println(f, join(vcat(f, ["instance", "nVertices", "nEdges"], [[alg_labels; alg .* ref_labels] for (alg, alg_labels) in algorithms]...), ", "))
    for (example_name, G) in examples
      l = "$example_name, $(n_vertices(G)), $(n_edges(G)), "
      K = uniform_hypergraph(G)
      print(l)
      timings = []
      for (algo, labels) in algorithms
        result = run_function(run_benchmark, K, algo, QQ; remote=useremote)
        # If the process died for some reason, put in appropriate number of "oom" or "oot".
        # The function might also have put n/a in the ref! columns; cf the las vegas algorithm.
        n_columns = length(labels) + length(ref_labels)
        append!(timings, isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result)
      end
      l *= join(timings, ", ")
      println(f, l)
      flush(f)
      println(l)
    end
  end
end

@info "Benchmarking Bipartite Graph examples"
bipartite_table(examples)

algorithms2 = [algorithms[3], algorithms[6]] # only lv-algorithm

function surfaces_table()
  example_dir = joinpath(root_of_project, "examples")
  surfaces_dir = joinpath(example_dir, "surfaces")
  fields = [
    QQ => "QQ",
    fpField(UInt(2)) => "Ftwo",
    fpField(UInt(3)) => "Fthree",
    fpField(UInt(5)) => "Ffive",
    fpField(UInt(7919)) => "Ftgig",
  ]
  open(surfaces_table_path, "w") do f
    println(f, join(["instance", "dim", "nVertices", "nFaces", "orientable", "genus", "index", "q", ["$(f[2])_$(l)" for f in fields, (algo, label) in algorithms2 for l in add_ref_labels(algo, label)]...], ", "))
    for example_file in readdir(surfaces_dir)
      println("Surface: $example_file")
      K = load(joinpath(surfaces_dir, example_file))
      # Parse the filename into the parameters
      dim, n, orientable, genus, index = match(r"d(\d+)_n(\d+)_o(\d+)_g(\d+)_#(\d+)", example_file).captures
      # Process the dimensions of the complex separately
      for q in 0:tryparse(Int, dim)
        S = uniform_hypergraph(K, q+1)
        timings = [example_file, dim, n_vertices(S), length(faces(S)), orientable, genus, index, q]
        for (F, _) in fields, (algo, labels) in algorithms2
          result = run_function(run_benchmark, S, algo, F; remote=true, time_limit=0.5)
          n_columns = length(labels) + length(ref_labels)
          append!(timings, isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result)
        end
        println(f, join(timings, ", "))
        flush(f)
        println(join(timings, ", "))
      end
    end
  end
end

@info "Benchmarking Surfaces examples"
surfaces_table()

function non_surfaces_table()
  example_dir = joinpath(root_of_project, "examples")
  non_surfaces_dir = joinpath(example_dir, "non_surfaces")
  fields = [
    QQ => "QQ",
    fpField(UInt(2)) => "Ftwo",
    fpField(UInt(3)) => "Fthree",
    fpField(UInt(5)) => "Ffive",
    fpField(UInt(7919)) => "Ftgig",
  ]
  open(non_surfaces_table_path, "w") do f
    println(f, join(["instance", "nVertices", "nFaces", "H_1",
                     ["$(f[2])_$(l)" for f in fields, (algo, label) in algorithms2 for l in add_ref_labels(algo, label)]...], ", "))
    for example_file in readdir(non_surfaces_dir)
      println("Non Surface: $example_file")
      K = load(joinpath(non_surfaces_dir, example_file))
      for q in 1:dim(K)
        S = uniform_hypergraph(K, q+1)
        timings = [example_file, n_vertices(S), length(faces(S)), homology(K, 1)]
        for (F, _) in fields, (algo, labels) in algorithms2
          result = run_function(run_benchmark, S, algo, F; remote=true, time_limit=3)
          # Append result to timings; if process died for some reson, put in appropriate number of "n/a"
          append!(timings, isnothing(result) ? fill("oom", length(labels)+length(ref_labels)) : result)
        end
        println(f, join(timings, ", "))
        flush(f)
        println(join(timings, ", "))
      end
    end
  end
end
@info "Benchmarking non Surfaces examples"
non_surfaces_table()

# clean up
map(rmprocs, workers())
