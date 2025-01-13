using Oscar, DataStructures, Distributed

add_ref_labels(alg, alg_labels) = [alg_labels; alg .* ref_labels]

workers_ready = Queue{Task}()

# Asynchronously start a new worker process and initialize it.
function initialize_new_worker()
  global workers_ready
  future = @async begin
    pids = addprocs(1)
    @everywhere pids include("benchmark_helpers.jl")
    @everywhere pid begin # precompile some stuff we might need
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

# Run the function `f` with `args...` on the next worker, and shutdown the process afterwards.
# time_limit is in hours
function run_function(f, args...; remote=true, time_limit=1)
  if remote
    initialize_new_worker() # Initialize new worker asynchronously for later use
    pid = get_worker() # Get initialized worker to run f
    try
      future = @async remotecall_fetch(f, pid, args...) # call remotely on worker
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
function run_benchmark(K, algorithm, F; finite_field_lv_trials=100)
  # run(`prlimit -v$(10e9) -c0 --pid $(getpid())`) # limit memory usage, limit in bytes
  Oscar.randseed!(1)
  n = n_vertices(K)

  # A inject logging to the ref! functions
  logger = Logger()
  logging_rref_cf(m) = ref_ff_rc_wrapper!(m; logger=logger)
  logging_rref_fl(m) = rref_lazy_pivots!(m; logger=logger, result_is_shifted_n_homogeneous=K.k-1)

  # Just to force compilation
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl)

  # The lv algorithm might not run ref! at all.
  logger[:ref] = fill("n/a", length(ref_labels))
  if algorithm == "av"
    R, x = polynomial_ring(K, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_cf)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "avf"
    R, x = polynomial_ring(K, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_fl)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "hv"
    p = perm(reverse(1:n))
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_cf)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "hvf"
    p = perm(reverse(1:n))
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_fl)
    return (t.time, t.bytes, logger[:ref]...)
  elseif algorithm == "lv"
    p = perm(reverse(1:n))
    trials = F isa QQField ? 1 : finite_field_lv_trials
    t = @timed exterior_shift(F, K, p; las_vegas_trials=trials, timed=true, (ref!)=logging_rref_cf)
    return (t.time, t.bytes, t.value[2]..., logger[:ref]...)
  elseif algorithm == "lvf"
    p = perm(reverse(1:n))
    trials = F isa QQField ? 1 : finite_field_lv_trials
    t = @timed exterior_shift(F, K, p; las_vegas_trials=trials, timed=true, (ref!)=logging_rref_fl)
    return (t.time, t.bytes, t.value[2]..., logger[:ref]...)
  else
    error("Unknown algorithm type")
  end
end
