import StaticArrays: SVector, MVector
import LinearAlgebra: rdiv!
import Base: front

export StencilRecurrence, StencilRecurrencePlan

"""
    StencilRecurrence{N, S, Tbuffer}

# Properties
- `stencil::SVector{S, CartesianIndex{N}}`: The relative index of the stencil. Can contain `(0,0)` (see `coef`)
- `coef::SVector{S, Function}`: The coefficient associated with each relative index. The one associated with `CartesianIndex(0,0)` refers to a constant added to that entry. The functions should take the index of the stencil and the eltype as input. Coefficients should be as accurate as possible, i.e. return either `Rational` or `BigFloat`.
- `buffer::Tbuffer`: a buffer to store temp results
- `offset::CartesianIndex{N}`: the top-left most stencil index where the current step starts at.
- `lastind::Int`: the recurrence terminates if `offset[end] > lastind`.
"""
mutable struct StencilRecurrence{N, S, Tbuffer<:AbstractArray} <: AbstractLinearRecurrence{N}
    stencil::SVector{S, CartesianIndex{N}}
    coef::SVector{S, Function}
    buffer::Tbuffer
    offset::MVector{N, Int} # marks the current step of the recurrence.
    lastind::Int
end
buffer(R::StencilRecurrence) = R.buffer

"""
    StencilRecurrencePlan{N, S, INIT} <: AbstractLinearRecurrencePlan

# Properties
- `stencil::SVector{S, CartesianIndex{N}}`: The relative index of the stencil. Can contain `(0,0)` (see `coef`)
- `coef::SVector{S, Function}`: The coefficient associated with each relative index. The one associated with `CartesianIndex(0,0)` refers to a constant added to that entry. The functions should take the index of the stencil and the eltype as input. Coefficients should be as accurate as possible, i.e. return either `Rational` or `BigFloat`.
- `init::INIT`: the function used for initial values.
- `size::Dims{N}`: the size of the whole array.
- `offset::SVector{N, Int}`: where the recurrence starts at.
"""
struct StencilRecurrencePlan{N, S, INIT} <: AbstractLinearRecurrencePlan
    stencil::SVector{S, CartesianIndex{N}}
    coef::SVector{S, Function}
    init::INIT
    size::Dims{N}
    offset::SVector{N, Int}
end
#StencilRecurrencePlan(stencil::SVector{S, CartesianIndex{N}}, coef::SVector{S, Function}, init, size::Dims{N}, offset::CartesianIndex{N}) where {N,S} = StencilRecurrencePlan{N, S, typeof(init)}(stencil, coef, init, size, offset) 
StencilRecurrencePlan(stencil, coef, init, size) = StencilRecurrencePlan(stencil, coef, init, size, -minimum(stencil))

function rdiv!(R::StencilRecurrence, x)
    rdiv!(R.buffer, x)
    R
end

function init(P::StencilRecurrencePlan; T=Float64, init=:default)
    if init == :default
        buffer = P.init(T, P.size)
    elseif init == :rand
        buffer = CircularArray(rand(complex(T), front(P.size)..., P.offset[end]+1))
    end
    StencilRecurrence(P.stencil, P.coef, buffer, MVector(P.offset), last(P.size)), eachslice(view(buffer.data, fill(:, ndims(buffer)-1)..., axes(buffer)[end][1:end-1]), dims=ndims(buffer))
end

function step!(R::StencilRecurrence{N}) where N
    slice = last(R.offset)
    if slice > R.lastind
        return nothing
    end
    stencil, coef, buffer = R.stencil, R.coef, R.buffer
    ind = CartesianIndex(R.offset...):CartesianIndex(last.(slicesupport(buffer, slice, dims=ndims(buffer)))..., slice)
    for i in ind
        v = zero(eltype(R.buffer))
        for (c, j) in zip(coef, stencil)
            if iszero(j)
                @inbounds v += c(Tuple(i)...)
            else
                k = i + j
                if checkbounds(Bool, buffer, front(k.I)..., 1)
                    @inbounds v += c(Tuple(i)...)*buffer[k]
                end
            end
        end
        buffer[i] = v
    end
    R.offset[end] += 1
    if N==1
        view(buffer, ind[1])
    else
        view(buffer, ind)
    end
end