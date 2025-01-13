include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
bipartite_table_path = joinpath(root_of_project, "bipartite_table_test_$computer_name.csv")

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

# If `useremote` is set, worker processes are prepared for remote execution asynchronously while the benchmark is running.
# To be safe from memory leaks etc., the workers are removed at the end of each benchmark.
# Start a few workers to have them ready for the benchmark.
useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

function bipartite_table(examples, algorithms)
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

# Algorithms to run, and the labels for the columns for that algorithm.
algorithms = [
  "av"   => ["avTime", "avMemory"],
  "hv"   => ["hvTime", "hvMemory"],
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  "avf"  => ["avfTime", "avfMemory"],
  "hvf"  => ["hvfTime", "hvfMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]

bipartite_table(examples, algorithms)

# clean up
map(rmprocs, workers())
