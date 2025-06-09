@testset "Moments" begin
    # Test with standard normal distribution
    f(x) = exp(-x^2 / 2)
    d = NumericallyIntegrable(f, (-Inf, Inf)) # Infinite support by default

    # Test moments against known values for standard normal
    @test isapprox(mean(d), 0.0, atol = 1e-7)
    @test isapprox(var(d), 1.0, atol = 1e-7)
    @test isapprox(std(d), 1.0, atol = 1e-7)
    @test isapprox(skewness(d), 0.0, atol = 1e-7)
    # Kurtosis is excess kurtosis (μ4/σ⁴ - 3)
    @test isapprox(kurtosis(d), 0.0, atol = 1e-7)

    # Test with a uniform distribution on [0, 2]
    f_uniform(x) = 1.0 * (0 <= x <= 2)
    d_uniform = NumericallyIntegrable(f_uniform, (0, 2))
    @test isapprox(mean(d_uniform), 1.0, atol = 1e-7)
    @test isapprox(var(d_uniform), 1 / 3, atol = 1e-7)
    @test isapprox(std(d_uniform), sqrt(1 / 3), atol = 1e-7)
    @test isapprox(skewness(d_uniform), 0.0, atol = 1e-7)
    @test isapprox(kurtosis(d_uniform), -1.2, atol = 1e-7) # Excess kurtosis for uniform is -6/5
end
