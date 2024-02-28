using PseudostableRecurrences
using Test

@testset "Examples" begin
    # Write your tests here.
    @testset "Linear recurrence sequence" begin
        stencil = (CartesianIndex(-3), CartesianIndex(-2), CartesianIndex(-1))
        coefs = (n -> 2//3, n -> -3, n -> 10//3)
        f_init(T) = [sqrt(T(2)), sqrt(T(2)), sqrt(T(2)), T(0)]
        P = StencilRecurrencePlan{Real}(stencil, coefs, f_init, (100,))
        @test stable_recurrence(P)[end] ≈ sqrt(2)
    end
end

@testset "Extensions" begin
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