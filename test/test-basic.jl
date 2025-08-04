using NumericalDistributions
using Distributions
using QuadGK
using Test

@testset "Basic functionality" begin
    # Test with a simple normal distribution truncated to [-2, 2]
    f(x) = exp(-x^2 / 2) * (abs(x) < 2)
    d = NumericallyIntegrable(f, (-2, 2))

    # is a distribution
    @test d isa ContinuousUnivariateDistribution
    @test NumericallyIntegrable(f, (-2, 2.0)) == NumericallyIntegrable(f, (-2.0, 2))

    # Test normalization
    total_prob = quadgk(x -> pdf(d, x), -2, 2)[1]
    @test total_prob ≈ 1.0

    # Test PDF properties
    @test pdf(d, -3) == 0.0  # outside support
    @test pdf(d, 3) == 0.0   # outside support
    @test pdf(d, 0) > pdf(d, 0.001) # peaks at zero
    @test pdf(d, 0) > 1 / sqrt(2π) # because of cutoff

    # Test CDF properties
    @test cdf(d, -3) == 0.0
    @test cdf(d, 3) == 1.0
    @test 0 < cdf(d, 0) < 1
    @test cdf(d, 1) > cdf(d, 0) # monotonicity
end

@testset "Integral function" begin
    # Test with a simple normal distribution truncated to [-2, 2]
    f(x) = exp(-x^2 / 2) * (abs(x) < 2)
    d = NumericallyIntegrable(f, (-2, 2))

    # Test that integral over entire support equals 1
    @test integral(d, -2, 2) ≈ 1.0

    # Test that integral over partial range matches CDF difference
    @test integral(d, 0.1, 0.5) ≈ cdf(d, 0.5) - cdf(d, 0.1)
    @test integral(d, -1, 1) ≈ cdf(d, 1) - cdf(d, -1)

    # Test boundary clamping - bounds outside support
    @test integral(d, -10, 10) ≈ 1.0  # should be clamped to [-2, 2]
    @test integral(d, -5, -3) ≈ 0.0   # entirely outside support
    @test integral(d, 3, 5) ≈ 0.0     # entirely outside support

    # Test partial boundary clamping
    @test integral(d, -10, 0) ≈ cdf(d, 0)    # lower bound clamped
    @test integral(d, 0, 10) ≈ 1.0 - cdf(d, 0)  # upper bound clamped

    # Test with different distribution - uniform on [0, 1]
    uniform_f(x) = (0 <= x <= 1) ? 1.0 : 0.0
    uniform_d = NumericallyIntegrable(uniform_f, (0, 1))
    
    @test integral(uniform_d, 0, 1) ≈ 1.0
    @test integral(uniform_d, 0.25, 0.75) ≈ 0.5
    @test integral(uniform_d, -1, 2) ≈ 1.0  # clamped to [0, 1]
end
