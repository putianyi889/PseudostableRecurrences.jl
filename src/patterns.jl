import Base: in, intersect, union, isempty
import Infinities: ℵ₀

abstract type EntryPattern end

struct EmptyPattern <: EntryPattern end
isempty(::EmptyPattern) = true
in(i, ::EmptyPattern) = false
intersect(::EmptyPattern, ::EntryPattern) = EmptyPattern()
intersect(::EntryPattern, ::EmptyPattern) = EmptyPattern()
intersect(::EmptyPattern, ::EmptyPattern) = EmptyPattern()

struct LatticePattern{EO,N} <: EntryPattern{N} end
LatticePattern(true) = LatticePattern{1}()
LatticePattern(false) = LatticePattern{0}()
in(i::CartesianIndex{N}, ::LatticePattern{1,N}) where N = isodd(sum(i.I))
in(i::CartesianIndex{N}, ::LatticePattern{0,N}) where N = iseven(sum(i.I))

abstract type AbstractBandedPattern <: EntryPattern{2} end
in(i::CartesianIndex{2}, P::AbstractBandedPattern) = upperbandwidth(P) ≤ i[1]-i[2] ≤ lowerbandwidth(P)
intersect(P::AbstractBandedPattern...) = BandedPattern(minimum(lowerbandwidth, P), minimum(upperbandwidth, P))

struct BandedPattern{T<:Integer} <: AbstractBandedPattern
    l::T
    u::T
end
upperbandwidth(P::BandedPattern) = P.u
lowerbandwidth(P::BandedPattern) = P.l

struct LowerBandedPattern{T<:Integer} <: AbstractBandedPattern
    l::T
end
upperbandwidth(::LowerBandedPattern) = ℵ₀
lowerbandwidth(P::LowerBandedPattern) = P.l
intersect(P::LowerBandedPattern...) = LowerBandedPattern(minimum(lowerbandwidth, P))

struct UpperBandedPattern{T<:Integer} <: AbstractBandedPattern
    u::T
end
upperbandwidth(P::UpperBandedPattern) = P.u
lowerbandwidth(::UpperBandedPattern) = ℵ₀
intersect(P::UpperBandedPattern...) = UpperBandedPattern(minimum(upperbandwidth, P))

struct StaticBandedPattern{L,U} <: AbstractBandedPattern end
upperbandwidth(::StaticBandedPattern{L,U}) where {L,U} = U
lowerbandwidth(::StaticBandedPattern{L,U}) where {L,U} = L

const TridiagonalPattern = StaticBandedPattern{1,1}
const DiagonalPattern = StaticBandedPattern{0,0}
const StaticLowerBandedPattern{L} = StaticBandedPattern{L,ℵ₀}
const StaticUpperBandedPattern{U} = StaticBandedPattern{ℵ₀,U}
const UpperTriangularPattern = StaticLowerBandedPattern{0}
const LowerTriangularPattern = StaticUpperBandedPattern{0}
const UpperHessenbergPattern = StaticLowerBandedPattern{1}
const LowerHessenbergPattern = StaticUpperBandedPattern{1}

abstract type AbstractDensePattern{N} <: EntryPattern{N} end
struct DensePattern{N,T<:Integer} <: AbstractDensePattern{N}
    lowerbounds::NTuple{N,T}
    upperbounds::NTuple{N,T}
end
DensePattern(A::AbstractArray) = DensePattern(map(first, axes(A)), map(last, axes(A)))

struct MarginPattern{N,T<:Integer} <: AbstractDensePattern{N}
    startmargin::NTuple{N,T}
    endmargin::NTuple{N,T}
end
struct FullPattern{N} <: AbstractDensePattern end

struct IntersectionPattern{P1<:EntryPattern,P2<:EntryPattern} <: EntryPattern
    p1::P1
    p2::P2
end

pattern(A::AbstractArray{T,N}) where {T,N} = FullPattern{N}()

findindex(A::AbstractArray, I...) = pattern_getindex(pattern(A), A, I...)
findindex(P::EntryPattern, A::AbstractArray, I...)

findindex(::FullPattern{N}, A::AbstractArray{T,N}, I...) where {T,N} = I
findindex(P::IntersectionPattern, A::AbstractArray, I)