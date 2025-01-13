# include("imports.jl")
using Oscar, DataStructures

Maybe{T} = Union{T, Nothing}
Logger = Dict{Symbol, Any}

function ref_ff_rc_wrapper!(m::AbstractAlgebra.Generic.MatSpaceElem{T}; logger::Maybe{Logger}=nothing) where T <: MPolyRingElem
  if isnothing(logger)
    return Oscar.ModStdQt.ref_ff_rc!(m)
  else
    s = _stats(m)
    r = Oscar.ModStdQt.ref_ff_rc!(m)
    append!(s, _stats(m))
    logger[:ref] = s
    return r
  end
end

"""
  rank(xs::Vector{Int}) -> Int

Rank simplex `xs` (ascendingly sorted list of 0-based vertices) in the (colexicographic) combinatorial number system.

# examples
```jldoctest
julia> rank([1,4,7,8])
59
``
"""
function rank(xs::Vector{Int})
	@req all(xs .> 0) "Can only decode positive numbers."
	@req issorted(xs) "xs must be sorted" #TODO Remove for efficiency?
	sum(binomial(v - 1, i) for (i, v) in enumerate(xs)) + 1
end

"""
  unrank(s::Int, k::Int) -> Vector{Int}

Unrank a `k`-combination from (colexicographic) combinatorial number system.

# examples
```jldoctest
julia> unrank(59, 4)
4-element Vector{Int64}:
 1
 4
 7
 8
```
"""
function decode(s::Int, k::Int)
	@req s > 0 "Can only decode positive numbers."
	s -= 1
	xs = Vector{Int}(undef, k) # The decoded vertices
	i = k + 1 # Index of the decoded vertex in xs
	while i >= 1
		v = 1::Int # The next decoded vertex is the largest v with binom(v, i) <= s
		while binomial(v - 1, i) <= s
			v += 1
		end
		v -= 1
		xs[i] = v
		s -= binomial(v - 1, i) # Now s = {s[1], ..., s[i-1]}
		i -= 1
	end
	return xs
end

"""
  encode_lex(xs::Vector{Int}, n::Int) -> Int

Rank a `k`-combination of {1,...,n} in the lexicographic combinatorial number system.
"""
function rank_lex(xs::Vector{Int}, n::Int)
	@req all(xs .> 0) "Can only decode positive numbers."
	@req issorted(xs) "xs must be sorted" #TODO Remove for efficiency?
  k = length(xs)
	binomial(n, k) - sum(binomial(n-v, k+1-i) for (i, v) in enumerate(xs))
end

"""
  decode_lex(s::Int, n::Int, k::Int) -> Vector{Int}

Unrank a `k`-combination of {1,...,n} from the lexicographic combinatorial number system.
"""
function unrank_lex(s::Int, n::Int, k::Int)
	@req s > 0 "Can only decode positive numbers."
  s = binomial(n, k)-s
	xs = Vector{Int}(undef, k) # The decoded vertices
  for i in k:-1:1 # Index of the decoded vertex in xs
		v = 0::Int # The next decoded vertex is the largest v with binom(v, i) <= s
		while binomial(v+1, i) <= s
			v += 1
		end
		xs[k+1-i] = n-v
		s -= binomial(v, i) # Now s = {s[1], ..., s[i-1]}
	end
	return xs
end

raw"""
  up_shift_indices(j::Int, n::Int, k::Int) -> Vector{Int}

Return the indices (in the lexicographic order on ``\binom{n}{k}``) of the simplices
that are optained from the simplex ranked `j` by increasing one vertex by one.

# Examples
```jldoctest
julia> n,k=7 ,2;
julia> unrank_lex.(7,n,k)
3-element Vector{Int64}:
 1
 3
 5
julia> unrank_lex.(up_shift_indices(7,n,k),n,k)
3-element Vector{Vector{Int64}}:
 [1, 3, 6]
 [1, 4, 5]
 [2, 3, 5]
```
"""
function up_shift_indices(j::Int, n::Int, k::Int)
  xs = unrank_lex(j, n, k)
  return [
    j + binomial(n-xs[i], k+1-i) - binomial(n-xs[i]-1, k+1-i) # replace vertex x[i] by x[i]+1...
    for i in k:-1:1
    if (
         (i  < k && xs[i]+1 != xs[i+1]) # ...if x[i]+1 is not already a vertex...
      || (i == k && xs[i]+1 <= n)       # ...and x[i]+1 is not out of bounds
    )
  ]
end

function rref_lazy_pivots!(m::AbstractAlgebra.Generic.MatSpaceElem{T}; logger::Maybe{Logger}=nothing, n::Int = 0, k::Int = 0) where T <: MPolyRingElem
  v = identity_matrix(base_ring(m), size(m, 1))
  r = 0
  I = Int[]
  s = isnothing(logger) ? [] : _stats(m)

  # cleared_columns = fill(false, size(m, 2))
  n_cleared = 0
  for j = 1:size(m, 2)
    # Evaluate j-th column of m * v
    c = v[r+1:end,:] * m[:,j:j]
    m[:,j:j] = zero(m[:,j:j])
    # @assert !cleared_columns[j] || c == zero(c) "Column $j should be cleared, but is non-zero."
    if iszero(c)
      # If we know that the result is shifted, no up-shift column of column j can have a step
      if n>0 && k>0
        # println("$j is non-face --> $(up_shift_indices(j, n, k)) are non-faces")
        cols = up_shift_indices(j, n, k)
        m[:, cols] .= zero(base_ring(m))
        n_cleared += length(cols)
        # cleared_columns[up_shift_indices(j, n, k)] .= true
      end
      continue
    end
    m[r+1,j] = one(base_ring(m))
    push!(I, j)
    # Break if this is the last necessary column
    if r+1 == size(m, 1)
      r += 1
      break
    end
    # Find the shortest non-zero polynomial in c
    _, i = efindmin(length, c[:,1]; filter=!iszero)
    # Use that as pivot; move corresponding row into row r+1
    if i > 1
      swap_rows!(v, r+i, r+1)
      swap_rows!(c,   i,   1)
    end
    # Eliminate other entries of c
    for i in 2:size(m, 1)-r
      if iszero(c[i])
        continue
      end
      _, a, b = gcd_with_cofactors(c[i], c[1])
      v[r+i,:] = b*v[r+i:r+i,:] - a*v[r+1:r+1,:]
      v[r+i,:] = divexact(v[r+i:r+i,:], content(v[r+i:r+i,:]))
    end
    r += 1
  end
  if !isnothing(logger)
    append!(s, _stats(m))
    logger[:ref] = s
  end
  if n != 0
    println(stderr, "cleared: $n_cleared/$(size(m, 2))")
  end
  return r
end

_stats(m::AbstractAlgebra.Generic.MatSpaceElem{T}) where T <: MPolyRingElem = [maximum(length.(m)), maximum(total_degree.(m))]

ref_labels = ["maxLenBefore", "maxDegBefore", "maxLenAfter", "maxDegAfter"]
