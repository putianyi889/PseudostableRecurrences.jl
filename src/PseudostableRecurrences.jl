module PseudostableRecurrences

# Write your package code here.

include("Utils.jl")
include("ArrayLayoutsExt.jl")
include("BandedMatricesExt.jl")
include("AbstractRecurrences.jl")
include("StencilRecurrences.jl")

"""
    precision_shift(P::AbstractLinearRecurrencePlan)

Estimates `log2` of the amplification of `P` by performing a full recurrence based on random initial conditions.
"""
function precision_shift(P::AbstractLinearRecurrencePlan)
    test = recurrence_init(P; init=:rand)
    shift = 1
    testmax = 1.0
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

function stable_recurrence(P::AbstractLinearRecurrencePlan, ::Type{T}) where T
    envprec = precision(BigFloat)
    varprec = precision(T)
    setprecision(varprec + Int(ceil(precision_shift(P))))
    R, ini = init(P, T=BigFloat, init=:default)
    container = [precision_convert.(T, varprec, y) for y in ini]
    while true
        v = step!(R)
        if isnothing(v)
            break
        end
        push!(container, precision_convert.(T, varprec, v))
    end
    setprecision(envprec)
    container
end

end
