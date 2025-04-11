"""
    NumericallyIntegrable{F,S} <: ContinuousUnivariateDistribution

A continuous univariate distribution defined by an unnormalized probability density function.

# Fields
- `unnormalized_pdf::F`: The unnormalized probability density function
- `integral::Float64`: The normalization constant (integral of the PDF)
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
struct NumericallyIntegrable{F,T<:Real} <: ContinuousUnivariateDistribution
    unnormalized_pdf::F
    integral::Float64
    support::Tuple{T,T}
    n_sampling_bins::Int

    function NumericallyIntegrable(f, support = (-Inf, Inf), n_sampling_bins = 300)
        integral = quadgk(f, support...)[1]
        F = typeof(f)
        T = eltype(support)
        new{F,T}(f, integral, support, n_sampling_bins)
    end
end

function Distributions.pdf(d::NumericallyIntegrable, x::Real)
    return d.support[1] < x < d.support[2] ? d.unnormalized_pdf(x) / d.integral : zero(x)
end

function Distributions.cdf(d::NumericallyIntegrable, x::Real)
    x ≤ d.support[1] && return zero(x)
    x ≥ d.support[2] && return one(x)
    integral, _ = quadgk(x -> d.unnormalized_pdf(x), d.support[1], x)
    return integral / d.integral
end

# Required method implementations for Distributions.jl interface
Distributions.minimum(d::NumericallyIntegrable) = d.support[1]
Distributions.maximum(d::NumericallyIntegrable) = d.support[2]
