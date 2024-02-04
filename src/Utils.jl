precision_convert(::Type{BigFloat}, prec, src::Real) = BigFloat(src, precision=prec)
precision_convert(T, prec, src::Real) = T(src)


"""
    slicetype(T)

Get the type of a slice of `T`

# Examples
```jldoctest
julia> using PseudostableRecurrences: slicetype

julia> slicetype(Matrix{Float64})
$(VERSION < v"1.6" ? "Array{Float64,1}" : "Vector{Float64} (alias for Array{Float64, 1})")

julia> slicetype(UnitRange{Int})
Int64

julia> slicetype(BitArray{3})
$(VERSION < v"1.6" ? "BitArray{2}" : "BitMatrix (alias for BitArray{2})")
```
"""
slicetype(T) = T
slicetype(T::Type{<:AbstractArray{S,N}}) where {S,N} = Union{Base.return_types(getindex, (T, fill(Colon,N-1)..., Int))...}

import CircularArrays: CircularArray
iscircular(::Array) = false
iscircular(::CircularArray) = true
iscircular(A::AbstractArray) = iscircular(parent(A))

tocircular(A::AbstractArray) = iscircular(A) ? A : CircularArray(A)

# not to be confused with LinearAlgebra.BLAS.dotu
_dotu(x, y) = mapreduce(*, +, x, y)