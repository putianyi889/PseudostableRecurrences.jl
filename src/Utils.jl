precision_convert(::Type{BigFloat}, prec, src::Real) = BigFloat(src, precision=prec)
precision_convert(T, prec, src::Real) = T(src)


"""
    slicetype(T)

Get the type of a slice of `T`

# Examples
```jldoctest
julia> using PseudostableRecurrences: slicetype # hide

julia> slicetype(Matrix{Float64})
$(
    if VERSION < v"1.6"
        "Array{Float64, 1}"
    else
        "Vector{Float64} (alias for Array{Float64, 1})"
    end
)

julia> slicetype(UnitRange{Int})
Int64

julia> slicetype(BitArray{3})
BitMatrix (alias for BitArray{2})
```
"""
slicetype(T) = T
slicetype(T::Type{<:AbstractArray{S,N}}) where {S,N} = Union{Base.return_types(getindex, (T, fill(Colon,N-1)..., Int))...}

import CircularArrays: CircularArray
iscircular(::Array) = false
iscircular(::CircularArray) = true
iscircular(A::AbstractArray) = iscircular(parent(A))

tocircular(A::AbstractArray) = iscircular(A) ? A : CircularArray(A)