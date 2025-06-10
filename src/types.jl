"""
    NumericallyIntegrable{F,S} <: ContinuousUnivariateDistribution

A continuous univariate distribution defined by an unnormalized probability density function.

# Fields
- `unnormalized_pdf::F`: The unnormalized probability density function
- `integral::I`: The normalization constant (integral of the PDF)
- `support::S`: The support range of the distribution (default: (-Inf, Inf))
- `n_sampling_bins::Int`: Number of bins used for sampling approximation

# Examples
```julia
# Define a custom distribution (e.g., truncated normal)
f(x) = exp(-x^2/2) * (abs(x) < 2)
dist = NumericallyIntegrable(f, (-2, 2))

# Sample from the distribution
samples = rand(dist, 1000)
```
"""
struct NumericallyIntegrable{F,T<:Real,I} <: ContinuousUnivariateDistribution
    unnormalized_pdf::F
    integral::I
    support::Tuple{T,T}
    n_sampling_bins::Int

    function NumericallyIntegrable(f, support; n_sampling_bins = 300)
        _integral = integral(f, support[1], support[2])
        F = typeof(f)
        T = eltype(support)
        I = typeof(_integral)
        new{F,T,I}(f, _integral, support, n_sampling_bins)
    end
end

"""
    integral(f, a::Real, b::Real)

Default numerical integration using adaptive Gauss-Kronrod quadrature via QuadGK.jl.
This method can be extended for custom function types to provide specialized integration algorithms,
as demonstrated in the README's custom integration example.
"""
integral(f, a::Real, b::Real) = quadgk(f, a, b)[1]

"""
    pdf(d::NumericallyIntegrable, x::Real)

Specialized PDF implementation that automatically handles normalization and support boundaries.
Divides the unnormalized PDF by the precomputed normalization integral and returns zero outside
the distribution's support range.
"""
function Distributions.pdf(d::NumericallyIntegrable, x::Real)
    return d.support[1] <= x <= d.support[2] ? d.unnormalized_pdf(x) / d.integral : zero(x)
end

"""
    cdf(d::NumericallyIntegrable, x::Real)

Dynamic CDF implementation that computes the cumulative probability through numerical integration.
Unlike most Distributions.jl types that have analytical CDFs, this method performs on-demand
integration of the unnormalized PDF from the lower support bound to `x`, normalized by the
total integral of the distribution.
"""
function Distributions.cdf(d::NumericallyIntegrable, x::Real)
    x ≤ d.support[1] && return zero(x)
    x ≥ d.support[2] && return one(x)
    _integral = integral(x -> d.unnormalized_pdf(x), d.support[1], x)
    return _integral / d.integral
end

# Required method implementations for Distributions.jl interface
"""
    minimum(d::NumericallyIntegrable)

Direct access to the lower bound of the distribution's support from the stored tuple.
"""
Distributions.minimum(d::NumericallyIntegrable) = d.support[1]

"""
    maximum(d::NumericallyIntegrable)

Direct access to the upper bound of the distribution's support from the stored tuple.
"""
Distributions.maximum(d::NumericallyIntegrable) = d.support[2]
