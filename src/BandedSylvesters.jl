import ArrayLayouts: colsupport, rowsupport
import BandedMatrices: bandwidth
import Base: size

export BandedSylvesterRecurrencePlan, BandedSylvesterRecurrence

"""
    BandedSylvesterRecurrence{T, TA<:AbstractMatrix{T}, TB<:AbstractMatrix{T}, TC<:AbstractMatrix{T}, TX<:AbstractMatrix{T}} <: AbstractLinearRecurrence{slicetype(TX)}

The recurrence generated from the infinite Sylvester equation ``AX+XB+C=0``, assuming ``X`` has infinite number of columns. The upper bandwidth of ``A`` has to be finite, the lower bandwidth of ``B`` has to be positive and finite and the lower bounding band of ``B`` can't contain zero. The recurrence is basically a cross-shaped stencil recurrence, where the width of the stencil is determined by The total bandwidth of ``B`` and the height by that of ``A``. See also [`BandedSylvesterRecurrencePlan`](@ref).

!!! note
    The restriction to the bandwidths may not be optimal, but it's the boundary of our knowledge at the moment.

# Properties
- `A::TA, B::TB, C::TC`: the matrices ``A``, ``B`` and ``C``. It's recommended to use `BandedMatrices.jl` to boost performance. Their dimensions don't have to match. Just make sure that no `BoundsError` happens during the recurrence.
- `buffer::TX`: the buffer that stores temp results.
- `sliceind::Int`: the current column index.
- `lastind::Int`: the last column index to be computed.
"""
BandedSylvesterRecurrence

@static if VERSION < v"1.8"
    include("../fragments/BandedSylvesterRecurrence1.jl")    
else
    include("../fragments/BandedSylvesterRecurrence2.jl")
end

buffer(R::BandedSylvesterRecurrence) = R.buffer
function rdiv!(P::BandedSylvesterRecurrence, x)
    rdiv!(P.buffer, x)
    P
end

"""
    BandedSylvesterRecurrencePlan{FA<:Function, FB<:Function, FC<:Function, INIT<:Function} <: AbstractLinearRecurrencePlan

See also [`BandedSylvesterRecurrence`](@ref).

# Properties
- `fA::FA, fB::FB, fC::FC`: the functions that generate `A`, `B` and `C` for the `BandedSylvesterRecurrence`. The functions should have the form `f(eltype, size...)`.
- `init::INIT`: the function that generates the first few columns of ``X`` in order to start the recurrence. It should have the form `f(eltype, size...)`.
- `size::Dims{2}`: the size of ``X``.
- `bandB::NTuple{2,Int}`: the bandwidths of `B`. Although `B` can have infinite upper bandwidth, that will cause the stencil to have infinite width and hence in practice, the upper bandwidth of `B` is always limited by the finite width of ``X``.
- `firstind::Int`: the first column to be computed. The default value is `bandB[2]+1`.
"""
struct BandedSylvesterRecurrencePlan{FA<:Function, FB<:Function, FC<:Function, INIT<:Function} <: AbstractLinearRecurrencePlan
    fA::FA
    fB::FB
    fC::FC
    init::INIT
    size::Dims{2}
    bandB::NTuple{2,Int}
    firstind::Int
end
BandedSylvesterRecurrencePlan(fA, fB, fC, init, size, bandB) = BandedSylvesterRecurrencePlan(fA, fB, fC, init, size, bandB, bandB[2]+1)
size(P::BandedSylvesterRecurrencePlan) = P.size

function init(P::BandedSylvesterRecurrencePlan; T=Float64, init=:default)
    if init == :default
        buffer = tocircular(method_to_precision(P.init, (Int,))(T, size(P)[1]))
    elseif init == :rand
        buffer = CircularArray(rand(T, size(P)[1], +(P.bandB...)))
    end
    A = method_to_precision(P.fA,())(T)
    B = method_to_precision(P.fB,())(T)
    C = method_to_precision(P.fC,())(T)
    BandedSylvesterRecurrence(A, B, C, buffer, P.firstind, P.size[end]), eachslice(view(buffer.data, fill(:, ndims(buffer)-1)..., axes(buffer)[end][1:P.firstind-1]), dims=ndims(buffer))
end

function step!(R::BandedSylvesterRecurrence)
    # A[m,:]X[:,n] + X[m,:]B[:,n] + C[m,n] = 0
    if R.sliceind > R.lastind
        return nothing
    end
    A, B, C, X = R.A, R.B, R.C, R.buffer
    nx = R.sliceind
    n = nx - bandwidth(B,1)
    ind = axes(X, 1)
    Bs = colsupport(B, n)
    I2 = (Base.OneTo(nx-1)) ∩ Bs
    Bv = B[I2, n]
    for m in ind
        I1 = rowsupport(A, m) ∩ ind
        X[m, nx] = -(_dotu(view(A,m,I1), view(X,I1,n)) + _dotu(view(X,m,I2), Bv) + C[m,n]) / B[nx,n]
    end
    R.sliceind += 1
    return view(X, :, nx)
end

