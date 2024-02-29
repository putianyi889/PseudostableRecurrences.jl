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
    @testset "A definite integral" begin
        stencil = (CartesianIndex(-1,-1), CartesianIndex(0,-1), CartesianIndex(0,0))
        coef1(m,n) = (m-1)//(n-1)
        coef2(m,n) = (m+n-3)//(n-1)//im
        coef3(m,n) = 2//(n-1)*ifelse(isodd(m+n), im, -1)
        coef = (coef1, coef2, coef3)
        function f_init(T, m)
            A = zeros(Complex{T},m,2)
            A[1,1] = π
            for mm in 3:2:m
                A[mm,1] = -A[mm-2,1] + 4//(mm-2)
            end
            A
        end
        P = StencilRecurrencePlan{Complex}(stencil, coef, f_init, (101,101), (1,2))
        stable_recurrence(P)
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