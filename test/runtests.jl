using Test
using NumericalDistributions
using Distributions
using QuadGK
using Statistics

@testset "NumericalDistributions.jl" begin
    @testset "Basic functionality" begin
        # Test with a simple normal distribution truncated to [-2, 2]
        f(x) = exp(-x^2 / 2) * (abs(x) < 2)
        d = NumericallyIntegrable(f, (-2, 2))

        # Test normalization
        total_prob, err = quadgk(x -> pdf(d, x), -2, 2)
        @test isapprox(total_prob, 1.0, rtol = 1e-6)

        # Test PDF properties
        @test pdf(d, -3) == 0.0  # outside support
        @test pdf(d, 3) == 0.0   # outside support
        @test pdf(d, 0) > pdf(d, 1)  # peak at center

        # Test CDF properties
        @test cdf(d, -3) == 0.0
        @test cdf(d, 3) == 1.0
        @test 0 < cdf(d, 0) < 1
        @test cdf(d, 1) > cdf(d, 0)  # monotonicity
    end

    @testset "Sampling" begin
        # Test with uniform distribution on [0,1]
        f(x) = 1.0
        d = NumericallyIntegrable(f, (0, 1))
        samples = rand(d, 10000)

        # Basic statistical tests
        @test all(0 .<= samples .<= 1)
        @test 0.45 < mean(samples) < 0.55  # should be close to 0.5
        @test 0.07 < std(samples) < 0.35   # should be close to 1/sqrt(12)
    end

    @testset "Infinite support" begin
        # Test with standard normal distribution
        f(x) = exp(-x^2 / 2)
        d = NumericallyIntegrable(f)

        # Test normalization
        total_prob, err = quadgk(x -> pdf(d, x), -Inf, Inf)
        @test isapprox(total_prob, 1.0, rtol = 1e-6)

        # Test sampling
        samples = rand(d, 10000)
        @test -0.1 < mean(samples) < 0.1    # should be close to 0
        @test 0.9 < std(samples) < 1.1      # should be close to 1
    end
end