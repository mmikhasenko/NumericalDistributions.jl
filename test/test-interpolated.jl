using NumericalDistributions
using Interpolations
using ForwardDiff
using QuadGK
using Test

a = interpolated(x -> 4 - x^2, -2:1.0:1; degree = Constant())
# plot(x -> pdf(a, x), -2, 2, fill = 0, fillalpha = 0.3, label = "PDF")
# plot(x -> cdf(a, x), -2, 2, fill = 0, fillalpha = 0.3, label = "CDF")
@testset "Interpolated constant" begin
    @test pdf(a, 0.4) ≈ 0.4705882369683137
    @test pdf(a, -0.5) ≈ 0.4705882369683137
    #
    @test cdf(a, -2.0) == 0.0
    @test cdf(a, -1.9) == 0.0
    @test cdf(a, 2.0) == 1.0
    @test cdf(a, 0.0) ≈ 0.5882352962103922
    #
    @test cdf(a, -1.9) ≈ quadgk(x -> pdf(a, x), -2, -1.9)[1]
    @test isapprox(cdf(a, 0.95), quadgk(x -> pdf(a, x), -2, 0.95)[1], atol = 1e-8)
    #
    @test a.integral ≈ 8.5
end

b = interpolated(x -> x^2, -2:0.5:1; degree = Linear())
# plot(x -> pdf(b, x), -2, 2, fill = 0, fillalpha = 0.3, label = "PDF")
# plot(y -> cdf(b, y), 0.9, 1.0, fillalpha = 0.3, label = "CDF")
# plot!(y -> quadgk(x -> pdf(b, x), -2, y)[1], 0.9, 1.0, label = "Integral", color = :red)
@testset "Interpolated linear" begin
    @test pdf(b, 0.4) ≈ 0.064
    #
    @test cdf(b, -2.0) == 0.0
    @test cdf(b, 2.0) == 1.0
    @test cdf(b, 0.0) ≈ 0.88
    #
    @test cdf(b, -1.9) ≈ quadgk(x -> pdf(b, x), -2, -1.9)[1]
    @test isapprox(cdf(b, 0.95), quadgk(x -> pdf(b, x), -2, 0.95)[1], atol = 1e-6)
    #
    @test b.integral ≈ 3.125
end

@testset "typeof vs eltype" begin
    g(x, p) = p * x
    h(p) = NumericallyIntegrable(x -> g(x, p), (0.0, 1.0)).integral
    @test isapprox(ForwardDiff.derivative(h, 2.0), 0.5; atol = 1e-12)
end
