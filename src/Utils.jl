precision_convert(::Type{BigFloat}, prec, src::Real) = BigFloat(src, precision=prec)
precision_convert(::Type{T}, prec, src::Complex) where T = Complex(precision_convert(T,prec,real(src)), precision_convert(T,prec,imag(src)))
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

include("EltypeExtensions.jl")

"""
    ToPrecision{F}(f::F) <: Function

Create a function where the first argument specifies the returned [`precisiontype`](@ref):
    (f::ToPrecision)(T, args...) = to_precision(T, f.f(args...))

# Examples
```jldoctest; setup = :(import PseudostableRecurrences: ToPrecision)
julia> f() = π
f (generic function with 1 method)

julia> f(n) = 1//n + (n-1)//n*im
f (generic function with 2 methods)

julia> g = ToPrecision(f)
(::$(repr(ToPrecision)){typeof(f)}) (generic function with 1 method)

julia> g(Float64)
3.141592653589793

julia> g(BigFloat)
3.141592653589793238462643383279502884197169399375105820974944592307816406286198

julia> g(Float16, 7)
Float16(0.1428) + Float16(0.857)im

julia> g(Float32, 3)
0.33333334f0 + 0.6666667f0im
```
"""
struct ToPrecision{F} <: Function
    f::F
end
(f::ToPrecision)(T, args...) = precisionconvert(precisiontype(T), f.f(args...))

"""
    method_to_precision(f, argtypes)

If `hasmethod(f, argtypes)`, return [`ToPrecision`](@ref)`(f)`. Otherwise return `f`.

# Examples
julia> f() = π

julia> Tf = method_to_precision(BigFloat, f, Tuple{})
"""
method_to_precision(f, argtypes) = hasmethod(f, argtypes) ? ToPrecision(f) : f

struct MaybeType{T,F} <: Function
    f::F
end
MaybeType{T}(f::Function) where T = MaybeType{T,typeof(f)}(f)
(f::MaybeType)(args...) = f(args)
function (f::MaybeType{T})(args) where T
    if applicable(f.f, args..., T)
        f.f(args...,T)
    else
        f.f(args...)
    end
end