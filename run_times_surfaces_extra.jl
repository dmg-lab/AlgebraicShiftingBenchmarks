
include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
surfaces_extra_table_path = joinpath(root_of_project, "surfaces_extra_test_$computer_name.csv")

# If `useremote` is set, worker processes are prepared for remote execution asynchronously while the benchmark is running.
# To be safe from memory leaks etc., the workers are removed at the end of each benchmark.
# Start a few workers to have them ready for the benchmark.
useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

function surfaces_table_extra(algorithms)
  example_dir = joinpath(root_of_project, "examples")
  surfaces_dir = joinpath(example_dir, "surfaces")
  fields = [
    0 => "QQ",
    2 => "F2",
    # 4 => "F4",
    3 => "F3",
    # 9 => "F9",
    5 => "F5",
    # 25 => "F25",
    7919 => "F7919",
		# 62710561 => "F62710561"
  ]
  instance_files = readdir(surfaces_dir)
  time_limit = 1/6
  println(stderr, "Benchmark can take up to $(time_limit * length(instance_files) * length(fields) * length(algorithms)) hours")
  open(surfaces_extra_table_path, "w") do f
    println(f, join(["instance", "dim", "nVertices", "nFaces", "orientable", "genus", "index", "q", ["$(f[2])$(uppercasefirst(l))" for f in fields, (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
    for example_file in readdir(surfaces_dir)
      println("Surface: $example_file")
      K = load(joinpath(surfaces_dir, example_file))
      # Parse the filename into the parameters
      dim, n, orientable, genus, index = tryparse.(Int, match(r"d(\d+)_n(\d+)_o(\d+)_g(\d+)_#(\d+)", example_file).captures)
      if n != 8 || orientable != 1 || genus != 0 || index <= 6
        println(stderr, "Skipping")
        continue
      end
      # Process the dimensions of the complex separately
      for q in 2:dim
        S = uniform_hypergraph(K, q+1)
        timings = [example_file, dim, n_vertices(S), length(faces(S)), orientable, genus, index, q]
        for (fieldsize, _) in fields, (algo, labels) in algorithms
          result = run_function(run_benchmark, S, algo, fieldsize; remote=true, time_limit=time_limit)
					n_columns = length(labels) + length(ref_labels)
          t = isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result
					append!(timings, t)
        end
        println(f, join(timings, ", "))
        flush(f)
      end
    end
  end
end

@info "Benchmarking surfaces extra"

# Algorithms to run, and the labels for the columns for that algorithm.
algorithms = [
  "av"   => ["avTime", "avMemory"],
  "hv"   => ["hvTime", "hvMemory"],
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  "avf"  => ["avfTime", "avfMemory"],
  "hvf"  => ["hvfTime", "hvfMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]

surfaces_table_extra(algorithms)

# clean up
map(rmprocs, workers())
