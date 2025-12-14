# See run_times_bijective.jl for a similar example with comments.

include("benchmark_helpers.jl")

const computer_name = gethostname()
const root_of_project = dirname(Base.active_project())
surfaces_extra_table_path = joinpath(root_of_project, "surfaces_extra_test_$computer_name.csv")

algorithms = [
  #"hv"   => ["hvTime", "hvMemory"],
  "lv" => ["lvTime", "lvMemory", "lvTrials", "lvTimeA", "lvMemoryA", "lvTimeB", "lvMemory"],
  #"hvf"  => ["hvfTime", "hvfMemory"],
  "lvf"  => ["lvfTime", "lvfMemory", "lvfTrials", "lvfTimeA", "lvfMemoryA", "lvfTimeB", "lvfMemory"]
]
fields = [
  4 => "F4"
]

useremote = true
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()
useremote && initialize_new_worker()

@info "Benchmarking surfaces extra"

surfaces_dir = joinpath(root_of_project, "examples", "surfaces")
instance_files = readdir(surfaces_dir)
time_limit = 1

println(stderr, "Benchmark can take up to $(time_limit * length(instance_files) * length(fields) * length(algorithms)) hours")
open(surfaces_extra_table_path, "w") do f
  # Table headers
  println(f, join(["#", "instance", "dim", "nVertices", "nFaces", "orientable", "genus", "index", "q", ["$(field)$(uppercasefirst(l))" for (_,field) in fields for (algo, label) in algorithms for l in add_ref_labels(algo, label)]...], ", "))
  for example_file in readdir(surfaces_dir)[12:34]
    K = load(joinpath(surfaces_dir, example_file))
    # Parse the filename into the parameters
    nr, dim, n, orientable, genus, index = match(r"^(\d\d)_manifold_lex_d(\d)_n(\d)_o(\d)_g(\d)_(\d\d)\..*$", example_file).captures
    # Only compute the shifts of the surfaces that remain after the following:
    println("Surface: $example_file")
    # Produce timings for each field and algorithm
    prev_uhg = Dict{Tuple{Int, String}, UniformHypergraph}((F, algo.first) => uniform_hypergraph(Vector{Int}[]) for algo in algorithms for (F, _) in fields)
    for q in 1:tryparse(Int, dim)
      S = uniform_hypergraph(K, q+1)
      timings = [nr, example_file, dim, n_vertices(S), length(faces(S)), orientable, genus, index, q]
      for (fieldsize, _) in fields, (algo, labels) in algorithms
        result = run_function(run_benchmark, S, algo, fieldsize; remote=useremote, time_limit=time_limit, lower_uhg=prev_uhg[(fieldsize, algo)])
        # Append the results to the timings, or, if computation died or timed out, append correct number of "oom" or "oot" respectively.
        n_columns = length(labels) + length(ref_labels)
        if !isnothing(result) && !(result isa Symbol)
          value = popfirst!(result)
          if !isnothing(value)
            prev_uhg[(fieldsize, algo)] = value
          end
        end

        append!(timings, isnothing(result) ? fill("oom", n_columns) : result == :timed_out ? fill("oot", n_columns) : result)
      end
      println(f, join(timings, ", "))
      flush(f)
    end
  end
end

# clean up
map(rmprocs, workers())
