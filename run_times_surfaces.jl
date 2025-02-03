# See run_times_bijective.jl for a similar example with comments.

include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
surfaces_table_path = joinpath(root_of_project, "surfaces_table_test_$computer_name.csv")

algorithms = [
  "lv"   => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
]
fields = [
  0 => "QQ",
  2 => "F2",
  4 => "F4",
  3 => "F3",
  9 => "F9",
  5 => "F5",
  25 => "F25",
  7919 => "F7919",
  62710561 => "F62710561"
]

useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

@info "Benchmarking Surfaces examples"

surfaces_dir = joinpath(root_of_project, "examples", "surfaces")
instance_files = readdir(surfaces_dir)
time_limit = 0.5
println(stderr, "Benchmark can take up to $(time_limit * length(instance_files) * length(fields) * length(algorithms)) hours")
open(surfaces_table_path, "w") do f
  # Table headers
  println(f, join(["#, instance", "dim", "nVertices", "nFaces", "orientable", "genus", "index", "q", ["$(f[2])$(uppercasefirst(l))" for f in fields for (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
  for example_file in readdir(surfaces_dir)
    println("Surface: $example_file")
    K = load(joinpath(surfaces_dir, example_file))
    # Parse the filename into the parameters
    nr, dim, n, orientable, genus, index = match(r"^(\d\d)_manifold_lex_d(\d)_n(\d)_o(\d)_g(\d)_(\d\d)\..*$", example_file).captures
    # Only shift the hypergraph of the 2-facets of K
    for q in 2:tryparse(Int, dim)
      S = uniform_hypergraph(K, q+1)
      # Initial entries of the row
      timings = [nr, example_file, dim, n_vertices(S), length(faces(S)), orientable, genus, index, q]
      # Produce timings for each field and algorithm
      for (F, _) in fields, (algo, labels) in algorithms
        result = run_function(run_benchmark, S, algo, F; remote=useremote, time_limit=time_limit)
        # Append the results to the timings, or, if computation died or timed out, append correct number of "oom" or "oot" respectively.
        n_columns = length(labels) + length(ref_labels)
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


