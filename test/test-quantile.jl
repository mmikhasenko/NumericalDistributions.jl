using NumericalDistributions
using Interpolations
using Test

f(x) = (1.2 - x^2) * exp(-x / 2)
grid = range(-0.5, 1, length = 11)
# Example: uniform distribution on [0, 1]
d_const = interpolated(f, grid; degree = Constant())

@testset "Inverse CDF for InterpolatedConstant" begin
    # Test quantile(cdf(x)) ≈ x
    xs = range(0, 1, length = 61)[2:end-1]  # Avoid endpoints for constant interpolation
    for x in xs
        u = cdf(d_const, x)
        x2 = quantile(d_const, u)
        @test isapprox(x, x2; atol = 1e-6)
    end
    # Test quantile(0) and quantile(1)
    @test isapprox(quantile(d_const, 0), d_const.support[1]; atol = 1e-6)
    @test isapprox(quantile(d_const, 1), d_const.support[2]; atol = 1e-6)
end


d_linear = interpolated(f, grid; degree = Linear())

@testset "Inverse CDF for InterpolatedLinear" begin
    # Example: triangular distribution on [0, 1]
    xs = range(0, 1, length = 61)[2:end-1]  # Avoid endpoints for linear interpolation
    for x in xs
        u = cdf(d_linear, x)
        x2 = quantile(d_linear, u)
        @test isapprox(x, x2; atol = 1e-6)
    end
    @test isapprox(quantile(d_linear, 0), d_linear.support[1]; atol = 1e-6)
    @test isapprox(quantile(d_linear, 1), d_linear.support[2]; atol = 1e-6)
end
