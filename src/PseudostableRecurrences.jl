module PseudostableRecurrences

using EltypeExtensions
using EltypeExtensions: _to_precisiontype

include("Utils.jl")
include("CircularArraysExt.jl")
include("FillArraysExt.jl")
include("BandedMatricesExt.jl")
include("AbstractRecurrences.jl")
include("StencilRecurrences.jl")
include("BandedSylvesters.jl")

import LinearAlgebra: norm
export stable_recurrence, precision_shift

"""
    precision_shift(P::AbstractLinearRecurrencePlan)

Estimates `log2` of the amplification of `P` by performing a full recurrence based on random initial conditions.
"""
function precision_shift(P::AbstractLinearRecurrencePlan)
    test, start = init(P; init=:rand)
    shift = 4 # safe choice from experiments
    testmax = norm(start, Inf)
    rdiv!(test, testmax)
    while true
        v = step!(test)
        if isnothing(v)
            break
        end
        testmax = norm(v, Inf)
        if testmax > 1e300
            shift += log2(testmax)
            rdiv!(test, testmax)
            testmax = 1.0
        end
    end
    shift += log2(testmax)
    return shift
end

function stable_recurrence(P::AbstractLinearRecurrencePlan, ::Type{T} = Float64) where T
    envprec = precision(BigFloat)
    varprec = precision(T)
    setprecision(varprec + Int(ceil(precision_shift(P))))
    R, ini = init(P, T=BigFloat, init=:default)
    container = [convert_precisiontype(T, y, varprec) for y in ini]
    while true
        v = step!(R)
        if isnothing(v)
            break
        end
        push!(container, convert_precisiontype(T, v, varprec))
    end
    setprecision(envprec)
    container
end

end
