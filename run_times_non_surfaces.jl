include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
surfaces_table_path = joinpath(root_of_project, "non_surfaces_table_test_$computer_name.csv")

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
    QQ => "QQ",
    fpField(UInt(2)) => "Ftwo",
    fpField(UInt(3)) => "Fthree",
    fpField(UInt(5)) => "Ffive",
    fpField(UInt(7919)) => "Ftgig",
  ]
  open(non_surfaces_table_path, "w") do f
    println(f, join(["instance", "nVertices", "nFaces", "H_1",
                     ["$(f[2])_$(l)" for f in fields, (algo, label) in algorithms2 for l in add_ref_labels(algo, label)]...], ", "))
    for example_file in readdir(non_surfaces_dir)[1:1]
      println("Non Surface: $example_file")
      K = load(joinpath(non_surfaces_dir, example_file))
      for q in 1:dim(K)
        S = uniform_hypergraph(K, q+1)
        timings = [example_file, n_vertices(S), length(faces(S)), homology(K, 1)]
        for (F, _) in fields, (algo, labels) in algorithms
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

algorithms = [
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]

@info "Benchmarking non Surfaces examples"
non_surfaces_table(algorithms)

# clean up
map(rmprocs, workers())
