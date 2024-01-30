"""
    smartdiv(x, y)

Returns `x // y` when supported and `x / y` otherwise.
"""
smartdiv(x, y) = x / y
smartdiv(x::Integer, y::Integer) = x // y
smartdiv(x::Complex, y::Real) = complex(smartdiv(real(x), y), smartdiv(imag(x), y))
smartdiv(x::Number, y::Complex) = smartdiv(x*conj(y), abs2(y))

"""
    basetype(T)

Iterates `eltype` on `T` until convergence. Returns the limit.
"""
@generated function basetype(::Type{T}) where T
    t = T
    while true
        el = eltype(T)
        if el == t
            return el
        end
        t
    end
end

@generated function basetype_convert(::Type{S}, ::Type{T}) where {S,T}
    if eltype(S) == S
        return T
    end
    return nameof(S){basetype_convert(eltype(S), T)}
end

precision_convert(::Type{BigFloat}, prec, src::Real) = BigFloat(src, precision=prec)
precision_convert(T, prec, src::Real) = T(src)

slicetype(T) = eltype(T)
slicetype(T::Type{<:AbstractVector}) = eltype(T)
slicetype(T::Type{<:AbstractArray{S,N}}) where {S,N} = nameof(T){S, N-1}