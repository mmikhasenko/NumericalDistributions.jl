using NumericalDistributions
using Interpolations
using Distributions
using Statistics
using Random
using QuadGK
using Test
using ForwardDiff

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

@testset "BinnedDensity" begin
    # Test uniform distribution
    g(x) = 1.0
    lims = (0.0, 1.0)
    nBins = 100
    bd = BinnedDensity(g, lims, nBins)

    # Test structure
    @test length(bd.grid) == nBins + 1  # nBins + 1 points to define nBins intervals
    @test bd.grid[1] ≈ 0.0
    @test bd.grid[end] ≈ 1.0
    @test length(bd.cumarr) == nBins + 1  # includes 0 at start
    @test bd.cumarr[1] ≈ 0.0
    @test bd.cumarr[end] ≈ 1.0

    # Test sampling
    rng = MersenneTwister(123)  # for reproducibility
    samples = [rand(rng, bd) for _ ∈ 1:10000]

    # Statistical tests for uniform distribution
    @test all(0 .<= samples .<= 1)
    @test 0.45 < mean(samples) < 0.55  # should be close to 0.5
    @test 0.28 < std(samples) < 0.30   # should be close to 1/√12 ≈ 0.289

    # Test non-uniform distribution (triangular)
    g_tri(x) = x
    bd_tri = BinnedDensity(g_tri, (0.0, 1.0), nBins)
    samples_tri = [rand(rng, bd_tri) for _ ∈ 1:10000]

    # Statistical tests for triangular distribution
    @test all(0 .<= samples_tri .<= 1)
    @test 0.63 < mean(samples_tri) < 0.70  # should be close to 2/3
    @test isapprox(median(samples_tri), 0.707, atol = 0.02)  # should be close to √(0.5)
end

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

a = Interpolated(x -> abs(x), -2:0.5:1; degree = Constant())
# plot(x -> pdf(a, x), -2, 2, fill = 0, fillalpha = 0.3, label = "PDF")
# plot(x -> cdf(a, x), -2, 2, fill = 0, fillalpha = 0.3, label = "CDF")
@testset "Interpolated constant" begin
    @test pdf(a, 0.4) ≈ 0.2
    @test cdf(a, -2.0) == 0.0
    @test cdf(a, 2.0) == 1.0
    @test pdf(a, -0.5) ≈ 0.2
    @test cdf(a, 0.0) ≈ 0.8
    @test a.integral ≈ 2.5
end

b = Interpolated(x -> abs(x), -2:0.5:1; degree = Linear())
# plot(x -> pdf(b, x), -2, 2, fill = 0, fillalpha = 0.3, label = "PDF")
# plot(x -> cdf(b, x), -2, 2, fill = 0, fillalpha = 0.3, label = "CDF")
@testset "Interpolated linear" begin
    @test pdf(b, 0.4) ≈ 0.16
    @test cdf(b, -2.0) == 0.0
    @test cdf(b, 2.0) == 1.0
    @test cdf(b, 0.0) ≈ 0.8
    @test b.integral ≈ 2.5
end
@testset "typeof vs eltype" begin
    g(x, p) = p * x
    h(p) = NumericallyIntegrable(x -> g(x, p), (0.0, 1.0)).integral
    @test isapprox(ForwardDiff.derivative(h, 2.0), 0.5; atol = 1e-12)
end

@testset "Convolution" begin
    Δ = 0.005
    x1 = -3:Δ:5
    x2 = -4:Δ:5
    uni_support, σ = (-0.5, 3.5), 0.3
    #
    yv_pdf1 = pdf.(Uniform(uni_support...), x1)
    yv_pdf2 = pdf.(Normal(0, σ), x2)

    # Analytical convolution: Uniform + Normal
    conv_analytical(x) =
        (cdf(Normal(uni_support[1], σ), x) - cdf(Normal(uni_support[2], σ), x)) /
        (uni_support[2] - uni_support[1])

    # Vector API
    conv_interp = convolve_pdf(yv_pdf1, yv_pdf2; Δ = Δ, t0_1 = first(x1), t0_2 = first(x2))

    # Test normalization - integral of the any pdf should be 1
    @test cdf(conv_interp, conv_interp.support[2]) ≈ 1.0

    # Test shape
    @test length(conv_interp.unnormalized_pdf.knots[1]) ==
          length(yv_pdf1) + length(yv_pdf2) - 1

    # Test against analytical at a few points (use nearest index for fair comparison)
    for x in range(-1, 4, length = 5)
        v = pdf(conv_interp, x)
        u = conv_analytical(x)
        @show v, u
        @test isapprox(v, u; atol = 4e-4)
    end
end

@testset "Convolution with Distributions" begin

    # Uniform and Normal distributions
    d1 = truncated(Uniform(-0.5, 3.5), -1, 4.0)
    d2 = truncated(Normal(0, 0.3), -2, 2)

    # Analytical convolution: Uniform + Normal
    conv_analytical(x) =
        (cdf(Normal(-0.5, 0.3), x) - cdf(Normal(3.5, 0.3), x)) / (3.5 + 0.5)

    # Use the new convolve_pdfs API
    conv_dist = convolve_pdfs(d1, d2; gridsize = 1000)

    # Test normalization
    @test cdf(conv_dist, conv_dist.support[2]) ≈ 1.0

    # Test against analytical at a few points
    for x in range(-1, 4, length = 5)
        v = pdf(conv_dist, x)
        u = conv_analytical(x)
        @test isapprox(v, u; atol = 4e-4)
    end

    # Test error for infinite support
    d2′ = Normal(0, 0.3)
    @test_throws ErrorException convolve_pdfs(d1, d2′; gridsize = 1000)
end

