# will be a new package
"""
    elconvert(T, A)

Similar to `convert(T, A)`, but `T` refers to the eltype.

# Examples
```jldoctest; setup = :(using PseudostableRecurrences: elconvert)
julia> elconvert(Float64, 1:10)
1.0:1.0:10.0

julia> typeof(elconvert(Float64, rand(Int, 3, 3)))
$(repr("text/plain", Matrix{Float64}))
```
"""
elconvert(::Type{T}, A::AbstractArray) where T = AbstractArray{T}(A)
elconvert(::Type{T}, A::AbstractRange) where T = T(first(A)):T(step(A)):T(last(A))
elconvert(::Type{T}, A::AbstractSet) where T = AbstractSet{T}(A)

"""
    _to_eltype(T, S)

Convert type `S` to have the `eltype` of `T`.
"""
_to_eltype(::Type{T}, ::Type{Array{S,N}}) where {T,S,N} = Array{T,N}
_to_eltype(::Type{T}, ::Type{Set}) where T = Set{T}
_to_eltype(::Type{T}, ::Type{S}) where {T,S} = Base.return_types(elconvert, (Type{T}, S))

nutype(x) = nutype(typeof(x))
nutype(T::Type) = throw(MethodError(nutype, T))

"""
    basetype(T::Type)

Recursively apply `eltype` to `T` until convergence.

# Examples
```jldoctest; setup = :(using PseudostableRecurrences: basetype)
julia> basetype(Matrix{BitArray})
Bool

julia> basetype(Vector{Set{Complex{Float64}}})
$(repr("text/plain", Complex{Float64}))

julia> basetype([1:n for n in 1:10])
Int64
```
"""
basetype(x) = basetype(typeof(x))
basetype(::Type{T}) where T = eltype(T) == T ? T : basetype(eltype(T))

"""
    _to_basetype(T::Type, S::Type)

Convert type `S` to have the [`basetype`](@ref) of `T`.
"""
_to_basetype(::Type{T}, ::Type{S}) where {T,S} = eltype(S) == S ? T : _to_eltype(_to_basetype(T, eltype(S)), S)

"""
    baseconvert(T::Type, A)

Similar to `convert(T, A)`, but `T` refers to the [`basetype`](@ref).
"""
baseconvert(::Type{T}, A::S) where {T,S} = convert(_to_basetype(T,S), A)

"""
    precisiontype(T::Type)

Returns the type that decides the precision of `T`. The difference from [`basetype`](@ref) is that `precisiontype` unwraps composite basetypes such as `Complex` and that `precisiontype` is not generalised.

# Examples
```jldoctest; setup = :(using PseudostableRecurrences: precisiontype)
julia> precisiontype(Complex{Float32})
Float32

julia> precisiontype(Matrix{ComplexF64})
Float64
```
"""
precisiontype(x) = precisiontype(typeof(x))
precisiontype(::Type{T}) where T<:Real = T
precisiontype(::Type{Complex{T}}) where T = T
precisiontype(::Type{T}) where T = eltype(T) == T ? throw(MethodError(precisiontype, T)) : precisiontype(basetype(T))

"""
    _to_precisiontype(T::Type, S::Type)

Convert type `S` to have the [`precisiontype`](@ref) of `T`. An exception is that if `T<:Integer`, then `Rational` will also be unwrapped.

# Examples
```jldoctest; setup = :(using PseudostableRecurrences: _to_precisiontype)
julia> _to_precisiontype(Float64, Complex{Rational{Int}})
$(repr("text/plain", Complex{Float64}))

julia> _to_precisiontype(BigFloat, Matrix{Complex{Bool}})
$(repr("text/plain", Matrix{Complex{BigFloat}}))

julia> _to_precisiontype(Int, Complex{Rational{BigInt}})
Complex{Rational{Int64}}
```
"""
_to_precisiontype(::Type{T}, ::Type{Complex}) where T = Complex{T}
_to_precisiontype(::Type{T}, ::Type{Complex{S}}) where {T,S} = Complex{_to_precisiontype(T,S)}
_to_precisiontype(::Type{T}, ::Type{<:Rational}) where T<:Integer = Rational{T}
_to_precisiontype(::Type{T}, ::Type{S}) where {T,S} = eltype(S) == S ? T : _to_eltype(_to_precisiontype(T, eltype(S)), S)

"""
    precisionconvert(T::Type, A, prec=precision(T))

Convert `A` to have the [`precisiontype`](@ref) of `T`. If `T` has adjustable precision such as `BigFloat`, the precision can be specified by `prec`, otherwise `prec` takes no effect.

# Examples
```jldoctest; setup = :(using PseudostableRecurrences: precisionconvert)
julia> precisionconvert(BigFloat, 1//3+im, 128)
0.3333333333333333333333333333333333333338 + 1.0im

julia> precisionconvert(Float16, [[m/n for n in 1:3] for m in 1:3])
3-element $(repr(Vector{Vector{Float16}})):
 [1.0, 0.5, 0.3333]
 [2.0, 1.0, 0.6665]
 [3.0, 1.5, 1.0]
```
"""
precisionconvert(T,A) = precisionconvert(T,A,precision(T))
precisionconvert(::Type{T}, A::S, prec) where {T,S} = convert(_to_precisiontype(T,S), A)
precisionconvert(::Type{BigFloat}, x::Real, prec) = BigFloat(x, precision=prec)
precisionconvert(::Type{BigFloat}, x::Complex, prec) = Complex(BigFloat(real(x), precision=prec), BigFloat(imag(x), precision=prec))
precisionconvert(::Type{BigFloat}, A, prec) = precisionconvert.(BigFloat, A, precision=prec)