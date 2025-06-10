"""
    interpolated(f, grid; degree = Linear())

Factory function that creates a distribution from a function using grid-based interpolation.
Allows creation of distributions with either constant or linear interpolation between grid points,
optimizing performance for PDF evaluation and sampling while maintaining accurate representation
of the original function.

# Arguments
- `f`: The function to interpolate, representing an unnormalized probability density function.
  Can be any callable object that accepts numeric input and returns a non-negative value.
- `grid::AbstractVector`: Vector of points where the function will be evaluated for interpolation.
  These points should cover the desired support of the distribution and be sorted in ascending order.
- `degree = Linear()`: Interpolation degree, either `Linear()` or `Constant()` from the
  Interpolations.jl package. Linear provides smoother PDFs but Constant can be more efficient for sampling.

# Returns
- `NumericallyIntegrable`: A distribution object with the specified interpolation scheme,
  normalized automatically and ready for PDF evaluation, CDF computation, and sampling.

# Performance Notes
- For sampling efficiency, `Constant()` generally outperforms `Linear()`
- For PDF accuracy, `Linear()` typically gives better results
- The number of grid points directly affects both accuracy and performance
- For distributions with sharp features, use more grid points in those regions

# Example
````julia
using NumericalDistributions
using Interpolations: Constant, Linear
using Plots

# Create a bimodal distribution
f(x) = exp(-(x-1)^2/0.2) + 0.5*exp(-(x+2)^2/0.3)
grid = range(-4, 4, length=200)

# Create with linear interpolation (smoother PDF)
d_linear = interpolated(f, grid, degree=Linear())

# Create with constant interpolation (faster sampling)
d_const = interpolated(f, grid, degree=Constant())

# Generate samples
samples_const = rand(d_const, 10000)
samples_linear = rand(d_linear, 10000)

# Compare distributions
xv = range(-4, 4, length=500)
yv_const = pdf.(Ref(d_linear), xv)
yv_linear = pdf.(Ref(d_linear), xv)
````
"""
function interpolated(f, grid::AbstractVector; degree = Linear())
    itr = interpolate((grid,), f.(grid), Gridded(degree))
    NumericallyIntegrable(itr, extrema(grid); n_sampling_bins = length(grid) - 1)
end

"""
    InterpolatedLinear

Type alias for a NumericallyIntegrable distribution that uses linear interpolation.
"""
const InterpolatedLinear = NumericallyIntegrable{
    Interpolations.GriddedInterpolation{T,N,TC,Gridded{ST},K},
    IT,
    IntT,
} where {T,N,TC,ST<:Interpolations.Linear,K,IT<:Real,IntT<:Real}

"""
    InterpolatedConstant

Type alias for a NumericallyIntegrable distribution that uses constant interpolation.
"""
const InterpolatedConstant = NumericallyIntegrable{
    Interpolations.GriddedInterpolation{T,N,TC,Gridded{ST},K},
    IT,
    IntT,
} where {T,N,TC,ST<:Interpolations.Constant,K,IT<:Real,IntT<:Real}


"""
    cdf(d::InterpolatedConstant, x::Real)

Specialized CDF implementation for constant-interpolated distributions using analytical integration.
Takes advantage of the piecewise constant nature of the PDF to efficiently compute
the CDF without numerical integration, by using the precomputed analytical integral formulas
from interpolate-integral.jl.
"""
function Distributions.cdf(d::InterpolatedConstant, x::Real)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]

    # Handle out-of-bounds
    x <= grid[1] && return zero(x)
    x >= grid[end] && return one(x)

    return integral(itr, grid[1], x) / d.integral
end


"""
    cdf(d::InterpolatedLinear, x::Real)

Optimized CDF computation for linearly-interpolated distributions using analytical formulas.
Takes advantage of the piecewise linear nature of the PDF to efficiently compute the CDF
using closed-form integrals of linear segments rather than general numerical integration methods,
significantly improving performance.
"""
function Distributions.cdf(d::InterpolatedLinear, x::Real)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]

    x <= grid[1] && return zero(x)
    x >= grid[end] && return one(x)

    return integral(itr, grid[1], x) / d.integral
end
