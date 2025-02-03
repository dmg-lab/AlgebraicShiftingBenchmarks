using Oscar, DataStructures, Distributed

add_ref_labels(alg, alg_labels) = [alg_labels; alg .* ref_labels]

workers_ready = Queue{Task}()

# Asynchronously start a new worker process and initialize it.
function initialize_new_worker()
  global workers_ready
  future = @async begin
    pids = addprocs(1)
    @everywhere pids include("benchmark_helpers.jl")
    @everywhere pids begin # precompile some stuff we might need
      collect(GF(2))
    end
    return pids[1]
  end
  enqueue!(workers_ready, future)
end

# Get the oldest worker from workers_ready, and wait for its initialization to finish and return the pid.
function get_worker()
  global workers_ready
  return fetch(dequeue!(workers_ready))
end

function show_result(f)
  return (args...; kwargs...) -> begin
    result = f(args...; kwargs...)
    println("$f($args; $kwargs) --> $(result)")
		return result
  end
end

# Run the function `f` with `args...` on the next worker, and shutdown the process afterwards.
# time_limit is in hours
function run_function(f, args...; remote=true, time_limit=1, kwargs...)
  if remote
    initialize_new_worker() # Initialize new worker asynchronously for later use
    pid = get_worker() # Get initialized worker to run f
    try
      future = @async remotecall_fetch(show_result(f), pid, args...; kwargs...) # call remotely on worker
      println(future)
      if timedwait(()->istaskdone(future), time_limit * 60*60) == :timed_out
        @warn "Remote worker $pid timed out"
        return :timed_out
      else
        return fetch(future)
      end
    catch e
      if e isa InterruptException
        @warn "Worker $pid receied $e; stopp all workers"
        rmprocs(workers())
        return
      else
        @warn "Worker $pid failed with exception $e"
        showerror(stderr, e)
        return nothing
      end
    finally
      # @info "remove worker $pid"
      @async rmprocs(pid)
    end
  else
    return f(args...)
  end
end

include("reduction.jl")
function run_benchmark(K::UniformHypergraph, algorithm, fsize::Int; finite_field_lv_trials::Int64=500)
  Oscar.randseed!(1)
  n = n_vertices(K)
  p = perm(reverse(1:n))
  F = fsize == 0 ? QQ : (is_prime(fsize) ? fpField(UInt(fsize)) : GF(fsize))

  # Inject logging to the ref! functions
  logger = Logger()
  logging_rref_cf(m) = ref_ff_rc_wrapper!(m; logger=logger)
  logging_rref_fl(m) = rref_lazy_pivots!(m; logger=logger)

  # Just to force compilation
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl)

  # The lv algorithm might not run ref! at all.
  logger[:ref] = fill("n/a", length(ref_labels))

  # Run the respective algorithm
  if algorithm == "av"
    println("Running av algorithm")
    R, x = polynomial_ring(F, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_cf)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "avf"
    println("Running avf algorithm")
    R, x = polynomial_ring(F, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_fl)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "hv"
    println("Running hv algorithm")
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_cf)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "hvf"
    println("Running hvf algorithm")
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_fl)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "lv"
    println("Running lv algorithm")
    trials = (F isa QQField) ? 1 : finite_field_lv_trials
    t = @timed exterior_shift(F, K, p; las_vegas_trials=trials, timed=true, (ref!)=logging_rref_cf)
    return (t.time, t.bytes, t.value[2]..., logger[:ref]...)
  elseif algorithm == "lvf"
    println("Running lvf algorithm")
    trials = (F isa QQField) ? 1 : finite_field_lv_trials
    t = @timed exterior_shift(F, K, p; las_vegas_trials=trials, timed=true, (ref!)=logging_rref_fl)
    return (t.time, t.bytes, t.value[2]..., logger[:ref]...)
  else
    error("Unknown algorithm type")
  end
end
