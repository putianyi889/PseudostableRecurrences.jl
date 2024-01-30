import Base: getindex, size

struct Product{T,N,S<:NTuple{N,AbstractVector{T}}} <: AbstractArray{T,N}
    args::S
end
Product(A::AbstractVector...) = Product(A)

size(A::Product) = length.(A.args)
function getindex(A::Product, I::Integer...)
    @boundscheck checkbounds(A, I...)
    @inbounds getindex.(A.args, I)
end