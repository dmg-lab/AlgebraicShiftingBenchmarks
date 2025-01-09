# include("imports.jl")
using Oscar, DataStructures

Maybe{T} = Union{T, Nothing}
Logger = Dict{Symbol, Any}

function column_reduce!(m)
  F = base_ring(base_ring(m))
  pivot(column) = efindmin(length, column; filter=!iszero)[1]
  r = 0 # number of linearly independent columns discovered so far
  p = [0 for _ in 1:minimum(size(m))]
  for j in 1:size(m, 2)
    for i in 1:r-1
      m[i,j] == 0 && continue
      _, c, d = gcd_with_cofactors(m[i,j], m[i,p[i]])
      m[:,j] = d*m[:,j:j] - c*m[:,p[i]:p[i]]
      m[:,j] = divexact(m[:,j:j], content(m[:,j:j]))
    end
    i = efindmin(length, m[:,j]; filter=!iszero)[1]
    isnothing(i) && continue
    swap_rows!(m, i, j)
    (r += 1) >= size(m, 1) && break
    p[r] = j
  end
end


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

function rref_lazy_pivots!(m::AbstractAlgebra.Generic.MatSpaceElem{T}; logger::Maybe{Logger}=nothing) where T <: MPolyRingElem
  v = identity_matrix(base_ring(m), size(m, 1))
  r = 0
  I = Int[]
  s = isnothing(logger) ? [] : _stats(m)
  for j = 1:size(m, 2)
    # Evaluate j-th column of m * v
    c = v[r+1:end,:] * m[:,j:j]
    m[:,j:j] = zero(m[:,j:j])
    if iszero(c)
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
  return r
end

_stats(m::AbstractAlgebra.Generic.MatSpaceElem{T}) where T <: MPolyRingElem = [maximum(length.(m)), maximum(total_degree.(m))]

ref_labels = ["maxLenBefore", "maxDegBefore", "maxLenAfter", "maxDegAfter"]
