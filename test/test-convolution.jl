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
    conv_interp = fft_convolve(yv_pdf1, yv_pdf2; Δ = Δ, t0_1 = first(x1), t0_2 = first(x2))

    # Test normalization - integral of the any pdf should be 1
    @test cdf(conv_interp, conv_interp.support[2]) ≈ 1.0

    # Test shape
    @test length(conv_interp.unnormalized_pdf.knots[1]) ==
          length(yv_pdf1) + length(yv_pdf2) - 1

    # Test against analytical at a few points (use nearest index for fair comparison)
    for x in range(-1, 4, length = 5)
        v = pdf(conv_interp, x)
        u = conv_analytical(x)
        @test isapprox(v, u; atol = 4e-4)
    end
end

@testset "Convolution with Distributions" begin
    Δ = 0.01
    # Uniform and Normal distributions
    d1 = Uniform(-0.5, 3.5)
    d2 = truncated(Normal(0, 0.3), -2, 2)

    # Analytical convolution: Uniform + Normal
    conv_analytical(x) =
        (cdf(Normal(-0.5, 0.3), x) - cdf(Normal(3.5, 0.3), x)) / (3.5 + 0.5)

    # Use the new fft_convolve API
    conv_dist = fft_convolve(d1, d2; gridsize = 1000)

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
    @test_throws ErrorException fft_convolve(d1, d2′; gridsize = 1000)
end

@testset "generic_fft vs FFTW" begin
    using FFTW
    rng = MersenneTwister(22)
    for n in (1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 64)
        x = randn(rng, n)
        F_ref = FFTW.fft(x)
        F_gen = NumericalDistributions.generic_fft(x)
        @test isapprox(F_gen, F_ref; rtol = 1e-10, atol = 1e-10 * n)
        @test isapprox(
            NumericalDistributions.generic_ifft(F_gen),
            x;
            rtol = 1e-10,
            atol = 1e-10 * n,
        )
    end
end

@testset "Convolution backend equivalence" begin
    Δ = 0.005
    x1 = -3:Δ:5
    x2 = -4:Δ:5
    y1 = pdf.(Uniform(-0.5, 3.5), x1)
    y2 = pdf.(Normal(0, 0.3), x2)
    kw = (; Δ = Δ, t0_1 = first(x1), t0_2 = first(x2))
    _, h_fftw = NumericalDistributions._fft_convolve(y1, y2; kw...)
    _, h_gen = NumericalDistributions._generic_fft_convolve(y1, y2; kw...)
    _, h_dir = NumericalDistributions._direct_convolve(y1, y2; kw...)
    n = length(h_fftw)
    @test h_fftw ≈ h_gen rtol = 1e-9 atol = 1e-9 * n
    @test h_dir ≈ h_fftw
    @test h_dir ≈ h_gen
end

@testset "Small-grid convolution backends" begin
    Δ = 0.1
    x1 = 0:Δ:1
    x2 = 0:Δ:1
    y1 = pdf.(Uniform(0, 1), x1)
    y2 = pdf.(Uniform(0, 1), x2)
    kw_dir = (; Δ = Δ, t0_1 = first(x1), t0_2 = first(x2))
    kw_fft = (; kw_dir..., pow2 = false)
    _, h_fftw = NumericalDistributions._fft_convolve(y1, y2; kw_fft...)
    _, h_gen = NumericalDistributions._generic_fft_convolve(y1, y2; kw_fft...)
    _, h_dir = NumericalDistributions._direct_convolve(y1, y2; kw_dir...)
    @test h_fftw ≈ h_dir rtol = 1e-12
    @test h_gen ≈ h_dir rtol = 1e-12
end

@testset "Convolution with ReverseDiff" begin
    using ReverseDiff
    Δ = 0.01
    x1 = -1:Δ:1
    x2 = -1:Δ:1
    y1 = pdf.(Uniform(-0.5, 0.5), x1)
    σ = ReverseDiff.track(0.3)
    y2 = pdf.(Normal(0, σ), x2)
    kw = (; Δ = Δ, t0_1 = first(x1), t0_2 = first(x2))
    @test fft_convolve(y1, y2; kw...) isa NumericallyIntegrable
    @test fft_convolve(y1, y2; kw..., algorithm = :direct) isa NumericallyIntegrable
    @test_throws ErrorException fft_convolve(y1, y2; kw..., algorithm = :generic)
    _, h_auto = NumericalDistributions._convolve_vectors(y1, y2; kw...)
    _, h_dir = NumericalDistributions._convolve_vectors(y1, y2; kw..., algorithm = :direct)
    @test h_auto ≈ h_dir
    f(σ_vec) = begin
        σ_val = σ_vec[1]
        y2′ = pdf.(Normal(0, σ_val), x2)
        _, h = NumericalDistributions._convolve_vectors(y1, y2′; kw..., algorithm = :direct)
        sum(h)
    end
    σ0 = [0.3]
    @test isfinite(f(σ0))
    deriv = ReverseDiff.gradient(f, σ0)
    @test isfinite(deriv[1])
end

@testset "Convolution algorithm selection" begin
    y1 = rand(10)
    y2 = rand(8)
    @test NumericalDistributions._resolve_convolution_algorithm(
        y1,
        y2;
        algorithm = :auto,
        pow2 = true,
    ) == :fftw
    using ReverseDiff
    y2t = ReverseDiff.track.(y2)
    @test NumericalDistributions._resolve_convolution_algorithm(
        y1,
        y2t;
        algorithm = :auto,
        pow2 = true,
    ) == :direct
    @test_throws ErrorException fft_convolve(y1, y2t; Δ = 0.1, t0_1 = 0.0, t0_2 = 0.0, algorithm = :fftw)
end
