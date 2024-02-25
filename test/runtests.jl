using PseudostableRecurrences
using Test

@testset "PseudostableRecurrences.jl" begin
    # Write your tests here.
    @testset "Product" begin
        using PseudostableRecurrences: Product
        A = Product(1:2,1:3)
        @test A isa AbstractMatrix{Tuple{Int,Int}}
        @test A == [(m,n) for m in 1:2, n in 1:3]
    end
end

using Documenter
@testset "DocTest" begin
    doctest(PseudostableRecurrences; manual = true)
end

using Aqua
@testset "Project quality" begin
    Aqua.test_all(PseudostableRecurrences, ambiguities=false, piracies=true, deps_compat=true, unbound_args=false)
end