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
