using Oscar, DataStructures

Maybe{T} = Union{T, Nothing}
Logger = Dict{Symbol, Any}

# For a matrix with entries in a multivariate polynomial ring, return the maximal length and maximal degree of the polynomials.
_stats(m::AbstractAlgebra.Generic.MatSpaceElem{T}) where T <: MPolyRingElem = [maximum(length.(m)), maximum(total_degree.(m))]

# Column labels relevant for output tables.
ref_labels = ["maxLenBefore", "maxDegBefore", "maxLenAfter", "maxDegAfter"]

# Wrapper for the Gaussian row elimination in Oscar.ModStdQt.ref_ff_rc! that logs some statistics
# (maximal length, maximal degree) of the matrix before and after the elimination.
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

# Alternative (lazy) row reduction algorithm, as described in the paper.
function rref_lazy_pivots!(m::AbstractAlgebra.Generic.MatSpaceElem{T}; logger::Maybe{Logger}=nothing) where T <: MPolyRingElem
  v = identity_matrix(base_ring(m), size(m, 1))
  r = 0
  I = Int[]
  s = isnothing(logger) ? [] : _stats(m)

  for j = 1:size(m, 2)
    # Evaluate j-th column of m * v
    c = v[r+1:end,:] * m[:,j:j]
    m[:,j:j] = zero(m[:,j:j])
    iszero(c) && continue
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
  return r
end
