mutable struct BandedSylvesterRecurrence{T, TA<:AbstractMatrix{T}, TB<:AbstractMatrix{T}, TC<:AbstractMatrix{T}, TX<:AbstractMatrix{T}} <: AbstractLinearRecurrence{slicetype(TX)}
    A::TA
    B::TB
    C::TC
    buffer::TX
    sliceind::Int
    lastind::Int
end