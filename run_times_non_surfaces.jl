include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
non_surfaces_table_path = joinpath(root_of_project, "non_surfaces_table_test_$computer_name.csv")

# If `useremote` is set, worker processes are prepared for remote execution asynchronously while the benchmark is running.
# To be safe from memory leaks etc., the workers are removed at the end of each benchmark.
# Start a few workers to have them ready for the benchmark.
useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()


function non_surfaces_table(algorithms)
  example_dir = joinpath(root_of_project, "examples")
  non_surfaces_dir = joinpath(example_dir, "non_surfaces")
  fields = [
    0 => "QQ",
    2 => "F2",
    4 => "F4",
    3 => "F3",
    9 => "F9",
    5 => "F5",
    25 => "F25",
    7919 => "F7919"
  ]

  open(non_surfaces_table_path, "w") do f
    println(f, join(["instance", "nVertices", "nFaces", "Hone", ["$(f[2])$(uppercasefirst(l))" for f in fields, (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
    for example_file in readdir(non_surfaces_dir)
      println("Non Surface: $example_file")
      K = load(joinpath(non_surfaces_dir, example_file))
      for q in 1:dim(K)
        S = uniform_hypergraph(K, q+1)
        timings = [example_file, n_vertices(S), length(faces(S)), homology(K, 1)]
        for (F, _) in fields, (algo, labels) in algorithms
          result = run_function(run_benchmark, S, algo, F; remote=true, time_limit=3, finite_field_lv_trials=300)
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

algorithms = [
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  # "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]

@info "Benchmarking non Surfaces examples"
non_surfaces_table(algorithms)

# clean up
map(rmprocs, workers())
