"""
    AbstractRecurrencePlan

The abstract type of recurrence plans. A recurrence plan should include all the information that can perform a recurrence with a precision. It needs to support the following methods:

- `init(::AbstractRecurrencePlan; T=Float64, init=:default)`: generates an `AbstractRecurrence` object where the actual recurrence runs on. Returns a vector of the view of the initial steps as well.
  - `T` is the underlying type where the precision is implied. 
  - `init` tells how the initial values are generated. Must support `:default` for the forward recurrence to run. It is also suggested to support `:rand` for the pseudostablization algorithm.
"""
abstract type AbstractRecurrencePlan end

"""
    AbstractLinearRecurrencePlan <: AbstractRecurrencePlan
    
A recurrence plan where the results are linear w.r.t. the initial conditions. Pseudostablization algorithm on linear recurrences is implemented.
"""
abstract type AbstractLinearRecurrencePlan <: AbstractRecurrencePlan end

"""
    AbstractRecurrence{T}

The abstract type where a recurrence instance runs. It needs to support the following methods:

- `step!(::AbstractRecurrence)`: step forward the recurrence. Returns `nothing` if the recurrence terminates, otherwise returns a view of the stepping result.
"""
abstract type AbstractRecurrence{T} end

"""
    AbstractLinearRecurrence{T} <: AbstractRecurrence{T}

A recurrence where the results are linear w.r.t. the initial conditions. To support pseudostablization, the following methods are needed:

- `rdiv!(::AbstractLinearRecurrence, x::Number)`: divides the system by `x`. It is used to prevent floating point overflow.
- `LinearAlgebra.norm(step!(::AbstractLinearRecurrence), Inf)`: computes the ∞-norm of a step. It is used to get the amplification of the system.
"""
abstract type AbstractLinearRecurrence{T} <: AbstractRecurrence{T} end

