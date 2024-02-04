mutable struct BandedSylvesterRecurrence{T, TA<:AbstractMatrix{T}, TB<:AbstractMatrix{T}, TC<:AbstractMatrix{T}, TX<:AbstractMatrix{T}} <: AbstractLinearRecurrence{slicetype(TX)}
    const A::TA
    const B::TB
    const C::TC
    const buffer::TX
    sliceind::Int
    const lastind::Int
end