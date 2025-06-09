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
