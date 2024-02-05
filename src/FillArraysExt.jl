import Base: getindex, size

struct Product{S<:NTuple{N,AbstractVector}} <: AbstractArray{Tuple{eltype.(S)...},N}
    args::S
end
Product(A::AbstractVector...) = Product(A)

size(A::Product) = length.(A.args)
function getindex(A::Product, I::Integer...)
    @boundscheck checkbounds(A, I...)
    @inbounds getindex.(A.args, I)
end
function getindex(A::Product, I...)
    @boundscheck checkbounds(A, I...)
    @inbounds Product((getindex.(A.args, _tovecind.(I))))
end

_tovecind(i::AbstractVector) = i
_tovecind(i::Colon) = i
_tovecind(i::Integer) = [i]
_tovecind(i::Tuple) = i
_tovecind(i::CartesianIndex) = [i]