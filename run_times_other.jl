# See run_times_bijective.jl for a similar example with comments.

include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
other_examples_table_path = joinpath(root_of_project, "other_examples_table_test_$computer_name.csv")

algorithms = [
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
]
fields = [
  0 => "QQ",
  2 => "F2",
  4 => "F4",
  3 => "F3",
  9 => "F9",
#  5 => "F5",
#  25 => "F25",
#  7919 => "F7919",
#  62710561 => "F62710561"
]

useremote = false
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

@info "Benchmarking other examples"
other_examples_dir = joinpath(root_of_project, "examples", "other_examples")
open(other_examples_table_path, "w") do f
  # Table headers
  println(f, join(["instance", "dim", "nVertices", "q", "Hq-1", "nFaces", ["$(f[2])$(uppercasefirst(l))" for f in fields for (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
  for example_file in readdir(other_examples_dir)
    println("Other Examples: $example_file")
    K = load(joinpath(other_examples_dir, example_file))
    prev_uhg = Dict{Tuple{Int, String}, UniformHypergraph}((F, algo.first) => uniform_hypergraph(Vector{Int}[]) for algo in algorithms for (F, _) in fields)
    for q in 1:dim(K)
      S = uniform_hypergraph(K, q+1)
      # Initial entries of the row
      timings = [example_file, dim(K), n_vertices(S), q, homology(K, q-1), length(faces(S))]
      # Produce timings for each field and algorithm
      for (F, _) in fields, (algo, labels) in algorithms
        # result = run_function(run_benchmark, S, algo, F; remote=useremote, time_limit=3, finite_field_lv_trials=500, lower_uhg=prev_uhg[(F, algo)])
        # Append the results to the timings, or, if computation died or timed out, append correct number of "oom" or "oot" respectively.
        n_columns = length(labels) + length(ref_labels)
        result = [[uniform_hypergraph(Vector{Int}[], 0, 0)]; [:timed_out for _ in 1:n_columns]]
        if !isnothing(result) && !(result isa Symbol)
          value = popfirst!(result)
          if !isnothing(value)
            prev_uhg[(F, algo)] = value
          end
        end

        append!(timings, isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result)
      end
      println(f, join(timings, ", "))
      flush(f)
      println(join(timings, ", "))
    end
  end
end

# clean up
map(rmprocs, workers())
