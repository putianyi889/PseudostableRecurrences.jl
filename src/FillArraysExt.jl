import Base: getindex, size

struct Product{T,N,S<:Tuple} <: AbstractArray{T,N}
    args::S
end
Product(A::NTuple{N, AbstractVector}) where N = Product{Tuple{eltype.(fieldtypes(typeof(A)))...}, N, typeof(A)}(A)
Product(A::AbstractVector...) = Product(A)

size(A::Product) = length.(A.args)
function getindex(A::Product{T,N}, I::Vararg{Integer,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds getindex.(A.args, I)
end
function getindex(A::Product{T,N}, I::Vararg{Any, N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds Product((getindex.(A.args, _tovecind.(I))))
end

_tovecind(i::AbstractVector) = i
_tovecind(i::Colon) = i
_tovecind(i::Integer) = [i]
_tovecind(i::Tuple) = i
_tovecind(i::CartesianIndex) = [i]