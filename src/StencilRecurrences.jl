import StaticArrays: SVector, MVector
import LinearAlgebra: rdiv!
import Base: front, size
import LazyArrays: BroadcastArray

export StencilRecurrence, StencilRecurrencePlan

"""
    StencilRecurrence{N,T,S,
        COEF<:NTuple{S,AbstractArray{T,N}},
        TB<:AbstractArray{T,N},}
        (stencil, coef, buffer, slicestart, sliceend, lastind)

# Properties
For `coef` and `slicesupport`, tt's suggested to use lazy arrays for performance.
- `stencil::SVector{S, CartesianIndex{N}}`: The relative index of the stencil. Can contain `(0,0)` (see `coef`)
- `coef::COEF<:NTuple{S,AbstractArray{T,N}}`: The coefficient associated with each relative index. The one associated with `CartesianIndex(0,0)` refers to a constant added to that entry. It's suggested to use lazy arrays for performance.
- `buffer::CircularArray{T,N,TB<:AbstractArray{T,N}}`: a buffer to store temp results. 
- `slicesupport::TI<:AbstractVector{CartesianIndices{N}}`: a vector that contains the indices of unknown entries in each slice.
- `sliceind::Int`: the index of the current slice.
- `lastind::Int`: the recurrence terminates if `sliceind > lastind`.
"""
struct StencilRecurrence{N, T, S, COEF<:NTuple{S,AbstractArray{T}}, TB<:AbstractArray{T,N}} <: AbstractLinearRecurrence{slicetype(TB)}
    stencil::SVector{S, CartesianIndex{N}}
    coef::COEF
    buffer::CircularArray{T,N,TB}
    slicestart::MVector{N, Int}
    sliceend::MVector{N, Int}
    lastind::Int
end

"""
    StencilRecurrencePlan{N, S, COEF<:NTuple{S,Function}, INIT<:Function} <: AbstractLinearRecurrencePlan

# Properties
- `stencil::SVector{S, CartesianIndex{N}}`: The relative index of the stencil. Can contain `(0,0)` (see `coef`)
- `coef::COEF<:NTuple{S,Function}`: The coefficient associated with each relative index. The one associated with `CartesianIndex(0,0)` refers to a constant added to that entry. The functions should be in the form `f(I..., T)` where `I` is the index of the stencil and `T` is the suggested return type. Coefficients should be at least as accurate as `T`. Exact-value types such as `Irrational`, `Rational` or `Integer` would do the job, and if that's not possible, `BigFloat` would work as well.
- `init::INIT<:Function`: the function used for initial values. The functions should be in the form `f(I..., T)` where `I` is the size of the array and `T` is the eltype.
- `size::Dims{N}`: the size of the whole array.
- `offset::CartesianIndex{N}`: the very first index where the recurrence starts at.
"""
struct StencilRecurrencePlan{N, S, COEF<:NTuple{S,Function}, INIT<:Function} <: AbstractLinearRecurrencePlan
    stencil::SVector{S, CartesianIndex{N}}
    coef::COEF
    init::INIT
    size::Dims{N}
    offset::CartesianIndex{N}
end
#StencilRecurrencePlan(stencil::SVector{S, CartesianIndex{N}}, coef::SVector{S, Function}, init, size::Dims{N}, offset::CartesianIndex{N}) where {N,S} = StencilRecurrencePlan{N, S, typeof(init)}(stencil, coef, init, size, offset) 
StencilRecurrencePlan(stencil, coef, init, size) = StencilRecurrencePlan(stencil, coef, init, size, -minimum(stencil))

size(P::StencilRecurrencePlan) = P.size

function rdiv!(R::StencilRecurrence, x)
    rdiv!(R.buffer, x)
    R
end

function init(P::StencilRecurrencePlan; T=Float64, init=:default)
    if init == :default
        buffer = P.init(T, P.size)
    elseif init == :rand
        buffer = CircularArray(rand(T, front(P.size)..., P.offset[end]))
    end
    sliceind = P.offset[end]
    sliceend = MVector(last.(slicesupport(buffer, sliceind, dims=ndims(buffer)))..., sliceind)
    StencilRecurrence(P.stencil, (f->BroadcastArray{T}(splat(f), Product(axes(P)))).(P.coef), buffer, MVector(P.offset.I...), sliceend, P.size[end]), eachslice(view(buffer.data, fill(:, ndims(buffer)-1)..., axes(buffer)[end][1:end-1]), dims=ndims(buffer))
end

function step!(R::StencilRecurrence{N}) where N
    slice = last(R.slicestart)
    if slice > R.lastind
        return nothing
    end
    stencil, coef, buffer = R.stencil, R.coef, R.buffer
    ind = CartesianIndex(R.slicestart...):CartesianIndex(R.sliceend...)
    for i in ind
        v = zero(eltype(R.buffer))
        for (c, j) in zip(coef, stencil)
            if iszero(j)
                @inbounds v += c[Tuple(i)...]
            else
                k = i + j
                if checkbounds(Bool, buffer, front(k.I)..., 1)
                    @inbounds v += c[Tuple(i)...]*buffer[k]
                end
            end
        end
        buffer[i] = v
    end
    R.slicestart[end] += 1
    R.sliceend[end] += 1
    if N==1
        view(buffer, ind[1])
    else
        view(buffer, ind)
    end
end