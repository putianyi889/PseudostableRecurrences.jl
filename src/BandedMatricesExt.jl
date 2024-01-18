import BandedMatrices: AbstractBandedMatrix, bandwidths, inbands_getindex, inbands_setindex!
import Infinities: InfiniteCardinal, ℵ₀
import CircularArrays: CircularArray
import Base: setindex!, @propagate_inbounds, axes, size, checkbounds
import LinearAlgebra: rdiv!, ldiv!, rmul!, lmul!

# export PseudoBandedMatrix, LowerBandedMatrix, UpperBandedMatrix

struct PseudoBandedMatrix{T, D<:AbstractMatrix{T}, L<:Integer, U<:Integer} <: AbstractBandedMatrix{T}
    data::D
    l::L
    u::U
end
bandwidths(a::PseudoBandedMatrix) = (a.l, a.u)
axes(a::PseudoBandedMatrix) = axes(a.data)
size(a::PseudoBandedMatrix) = size(a.data)
checkbounds(a::PseudoBandedMatrix, I...) = checkbounds(a.data, I...)
@propagate_inbounds inbands_getindex(a::PseudoBandedMatrix, j::Integer, k::Integer) = a.data[j, k]
@propagate_inbounds inbands_setindex!(a::PseudoBandedMatrix, v, j::Integer, k::Integer) = setindex!(a.data, v, j, k)

const LowerBandedMatrix{T,D} = PseudoBandedMatrix{T,D,Int,InfiniteCardinal{0}}
const UpperBandedMatrix{T,D} = PseudoBandedMatrix{T,D,InfiniteCardinal{0},Int}
LowerBandedMatrix(data,l) = PseudoBandedMatrix(data,l,ℵ₀)
UpperBandedMatrix(data,u) = PseudoBandedMatrix(data,ℵ₀,u)

for op in (:rdiv!, :ldiv!, :rmul!, :lmul!)
    @eval function $op(a::PseudoBandedMatrix{<:Any, <:CircularArray}, x::Number)
        $op(a.data, x)
        a
    end
end