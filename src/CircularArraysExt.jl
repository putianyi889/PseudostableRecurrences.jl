struct SemiCircularArray{T,N,S<:AbstractArray{T,N},D} <: AbstractArray{T,N}
    data::S
    circulardims::NTuple{D,Int}
end