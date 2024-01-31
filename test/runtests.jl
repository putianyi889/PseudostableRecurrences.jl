using PseudostableRecurrences
using Test

@testset "PseudostableRecurrences.jl" begin
    # Write your tests here.
end

using Documenter
@testset "DocTest" begin
    doctest(PseudostableRecurrences; manual = true)
end

using Aqua
@testset "Project quality" begin
    Aqua.test_all(PseudostableRecurrences, ambiguities=false, piracies=true, deps_compat=false, unbound_args=false)
end