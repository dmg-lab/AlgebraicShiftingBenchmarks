# See run_times_bijective.jl for a similar example with comments.

include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
surfaces_extra_table_path = joinpath(root_of_project, "surfaces_extra_test_$computer_name.csv")

algorithms = [
  "hv"   => ["hvTime", "hvMemory"],
  # "lv" has been run already in run_times_surfaces.jl
  "hvf"  => ["hvfTime", "hvfMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]
fields = [
  0 => "QQ",
  2 => "F2",
  3 => "F3",
  5 => "F5",
  7919 => "F7919",
]

useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

@info "Benchmarking surfaces extra"

surfaces_dir = joinpath(root_of_project, "examples", "surfaces")
instance_files = readdir(surfaces_dir)
time_limit = 1/6
println(stderr, "Benchmark can take up to $(time_limit * length(instance_files) * length(fields) * length(algorithms)) hours")
open(surfaces_extra_table_path, "w") do f
  # Table headers
  println(f, join(["instance", "dim", "nVertices", "nFaces", "orientable", "genus", "index", "q", ["$(field)$(uppercasefirst(l))" for (_,field) in fields for (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
  for example_file in readdir(surfaces_dir)
    K = load(joinpath(surfaces_dir, example_file))
    # Parse the filename into the parameters
    dim, n, orientable, genus, index = tryparse.(Int, match(r"d(\d+)_n(\d+)_o(\d+)_g(\d+)_#(\d+)", example_file).captures)
    # Only compute the shifts of the surfaces that remain after the following:
    if n != 8 || orientable != 1 || genus != 0
      continue
    end
    println("Surface: $example_file")
      # Produce timings for each field and algorithm
    for q in 2:dim
      S = uniform_hypergraph(K, q+1)
      timings = [example_file, dim, n_vertices(S), length(faces(S)), orientable, genus, index, q]
      for (fieldsize, _) in fields, (algo, labels) in algorithms
        result = run_function(run_benchmark, S, algo, fieldsize; remote=true, time_limit=time_limit)
        # Append the results to the timings, or, if computation died or timed out, append correct number of "oom" or "oot" respectively.
        n_columns = length(labels) + length(ref_labels)
        append!(timings, isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result)
      end
      println(f, join(timings, ", "))
      flush(f)
    end
  end
end

# clean up
map(rmprocs, workers())
