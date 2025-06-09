using NumericalDistributions
using Interpolations
using Test
using QuadGK

itr_const = interpolate((1.0:3.0,), [1.0, 3.0, 2.0], Gridded(Constant()))

@testset "Interpolate Constant() integral" begin
    @test integral(itr_const, 1.0, 2.5) == 3.5
    # at the right edge of the grid
    @test integral(itr_const, 1.0, 3.0) == 4.5
    @test isapprox(
        integral(itr_const, 1.0, 2.999),
        integral(itr_const, 1.0, 3.0);
        atol = 1e-2,
    )
    # at then left edge of the grid
    @test integral(itr_const, 1.0, 1.1) ≈ 0.1
    # zero range
    @test integral(itr_const, 2.1, 2.1) == 0.0
    # integral over a small range
    @test integral(itr_const, 2.1, 2.2) ≈ 0.3
end

itr_linear = interpolate((1.0:3.0,), [1.0, 3.0, 2.0], Gridded(Linear()))

@testset "Interpolate Linear() integral" begin
    @test integral(itr_linear, 1.0, 2.5) ≈ quadgk(itr_linear, 1.0, 2.5)[1]
    # at the right edge of the grid
    @test integral(itr_linear, 1.0, 3.0) == quadgk(itr_linear, 1.0, 3.0)[1]
    @test isapprox(
        integral(itr_linear, 1.0, 2.999),
        integral(itr_linear, 1.0, 3.0);
        atol = 1e-2,
    )
    # at then left edge of the grid
    @test integral(itr_linear, 1.0, 1.1) ≈ quadgk(itr_linear, 1.0, 1.1)[1]
    # zero range
    @test integral(itr_linear, 2.1, 2.1) == 0.0
    # integral over a small range
    @test integral(itr_linear, 2.1, 2.2) ≈ quadgk(itr_linear, 2.1, 2.2)[1]
end

