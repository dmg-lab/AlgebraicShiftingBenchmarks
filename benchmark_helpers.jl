using Oscar, DataStructures, Distributed

function exterior_shift_lv_timed(F::Field, K::ComplexOrHypergraph, p::PermGroupElem; n_samples=100, kw...)
  # this might need to be changed based on the characteristic
  # we expect that the larger the characteristic the smaller the sample needs to be
  # setting to 100 now for good measure
  # Compute n_samples many shifts by radom matrices, and take the lexicographically minimal one, together with its first index of occurrence.
  
  random_matrices = [random_rothe_matrix(F, p) for _ in 1:n_samples]
  (shift, i), stats... = @timed Oscar.efindmin((exterior_shift(K, r) for (i, r) in enumerate(random_matrices)); lt=isless_lex)
  # Check if `shift` is the generic exterior shift of K
  prime_field = characteristic(F) == 0 ? QQ : fpField(UInt(characteristic(F)))
  n = n_vertices(K)
  # check shifted is algorithm 3 in the paper
  is_correct_shift, stats2... = @timed (p != perm(reverse(1:n)) || is_shifted(shift)) && Oscar.check_shifted(prime_field, K, shift, p; kw...)
  
  if is_correct_shift
    return shift, (i, stats.time, stats.bytes, stats2.time, stats2.bytes)
  else
    return nothing, (">$n_samples", stats.time, stats.bytes, stats2.time, stats2.bytes)
  end
end

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
    return f(args...; kwargs...)
  end
end

include("reduction.jl")
function run_benchmark(K::UniformHypergraph, algorithm, fsize::Int; finite_field_lv_trials::Int64=500, kw...)
  Oscar.randseed!(1)
  n = n_vertices(K)
  p = perm(reverse(1:n))
  F = fsize == 0 ? QQ : (is_prime(fsize) ? fpField(UInt(fsize)) : GF(fsize))

  # Inject logging to the ref! functions
  logger = Logger()
  logging_rref_cf(m) = ref_ff_rc_wrapper!(m; logger=logger)
  logging_rref_fl(m) = rref_lazy_pivots!(m; logger=logger)

  # Just to force compilation
  exterior_shift_lv_timed(QQ, uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_cf, las_vegas_trials=0)
  exterior_shift_lv_timed(QQ, uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl)
  exterior_shift(uniform_hypergraph([[1,3],[1,4]]); (ref!)=logging_rref_fl, las_vegas_trials=0)
  exterior_shift(klein_bottle())
  
  # The lv algorithm might not run ref! at all.
  logger[:ref] = fill("n/a", length(ref_labels))

  # Run the respective algorithm
  if algorithm == "av"
    println("Running av algorithm")
    R, x = polynomial_ring(F, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_cf, kw...)
    return [t.value, t.time, t.bytes, logger[:ref]...]
  elseif algorithm == "avf"
    println("Running avf algorithm")
    R, x = polynomial_ring(F, :x => (1:n, 1:n))
    g = matrix(R, x)
    t = @timed exterior_shift(K, g; (ref!)=logging_rref_fl, kw...)
    return [t.value, t.time, t.bytes, logger[:ref]...]
  elseif algorithm == "hv"
    println("Running hv algorithm")
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_cf, las_vegas_trials=0, kw...)
    return [t.value, t.time, t.bytes, logger[:ref]...]
  elseif algorithm == "hvf"
    println("Running hvf algorithm")
    t = @timed exterior_shift(F, K, p; (ref!)=logging_rref_fl, las_vegas_trials=0, kw...)
    return [t.value, t.time, t.bytes, logger[:ref]...]
  elseif algorithm == "lv"
    println("Running lv algorithm")
    trials = (F isa QQField) ? 1 : finite_field_lv_trials
    t = @timed exterior_shift_lv_timed(F, K, p; n_samples=trials, (ref!)=logging_rref_cf, kw...)
    return [t.value[1], t.time, t.bytes, t.value[2]..., logger[:ref]...]
  elseif algorithm == "lvf"
    println("Running lvf algorithm")
    trials = (F isa QQField) ? 1 : finite_field_lv_trials
    t = @timed exterior_shift_lv_timed(F, K, p; n_samples=trials, (ref!)=logging_rref_fl, kw...)
    return [t.value[1], t.time, t.bytes, t.value[2]..., logger[:ref]...]
  else
    error("Unknown algorithm type")
  end
end
