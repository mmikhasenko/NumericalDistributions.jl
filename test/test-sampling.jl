using NumericalDistributions.Interpolations
using NumericalDistributions
using Statistics
using Random
using QuadGK
using Test

@testset "Interpolated Sampling" begin
    # Test uniform distribution
    g(x) = 1.0
    lims = (0.0, 1.0)
    nBins = 100
    bd = interpolated(g, range(lims..., nBins + 1); degree = Constant())

    # Test sampling
    rng = MersenneTwister(123)  # for reproducibility
    samples = rand(rng, bd, 10000)

    # Statistical tests for uniform distribution
    @test all(0 .<= samples .<= 1)
    @test 0.45 < mean(samples) < 0.55  # should be close to 0.5
    @test 0.28 < std(samples) < 0.30   # should be close to 1/√12 ≈ 0.289

    # Test non-uniform distribution (triangular)
    g_tri(x) = x
    bd_tri = interpolated(g_tri, range(0.0, 1.0, nBins + 1); degree = Constant())
    samples_tri = rand(rng, bd_tri, 10000)

    # Statistical tests for triangular distribution
    @test all(0 .<= samples_tri .<= 1)
    @test 0.63 < mean(samples_tri) < 0.70  # should be close to 2/3
    @test isapprox(median(samples_tri), 0.707, atol = 0.02)  # should be close to √(0.5)
end

@testset "NumericallyIntegrable Sampling" begin
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
    d = NumericallyIntegrable(f, (-Inf, Inf))

    # Test normalization
    total_prob, err = quadgk(x -> pdf(d, x), -Inf, Inf)
    @test isapprox(total_prob, 1.0, rtol = 1e-6)

    # Test sampling
    samples = rand(d, 10000)
    @test -0.1 < mean(samples) < 0.1    # should be close to 0
    @test 0.9 < std(samples) < 1.1      # should be close to 1
end
