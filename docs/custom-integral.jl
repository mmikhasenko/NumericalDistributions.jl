using NumericalDistributions

struct FullPeriodSinCosSquared{T}
    which::T
end
(o::FullPeriodSinCosSquared)(x::Number) = o.which(x)^2  # here is the function to call
NumericalDistributions.integral(o::FullPeriodSinCosSquared, a::Real, b::Real) = (b - a) / 2

d_sin = NumericallyIntegrable(FullPeriodSinCosSquared(sin), (0, 2π))  # should work with customary integral
d_sin.integral # π

d_cos = NumericallyIntegrable(FullPeriodSinCosSquared(cos), (0, 2π))  # should work with customary integral
d_cos.integral # π

using Plots
theme(:boxed)
let
    plot()
    plot!(x -> pdf(d_sin, x), 0, 2π, label = "~sin²(x)", fill = 0, fillalpha = 0.3)
    plot!(x -> pdf(d_cos, x), 0, 2π, label = "~cos²(x)", fill = 0, fillalpha = 0.3)
end
