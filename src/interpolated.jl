"""
    interpolated(f, grid; degree = Linear())

Create a distribution using interpolation of a function over a grid of points.

# Arguments
- `f`: The function to be interpolated
- `grid::AbstractVector`: Grid points for interpolation
- `degree`: Interpolation type, either `Linear()` or `Constant()`

# Returns
- `NumericallyIntegrable`: A distribution based on the interpolated function
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

Compute the cumulative distribution function for a distribution with constant interpolation using an analytic expression for integral.

# Arguments
- `d::InterpolatedConstant`: The distribution
- `x::Real`: The point at which to evaluate the CDF

# Returns
- The CDF value at point x
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

Compute the cumulative distribution function for a distribution with linear interpolation using an analytic expression for integral.

# Arguments
- `d::InterpolatedLinear`: The distribution
- `x::Real`: The point at which to evaluate the CDF

# Returns
- The CDF value at point x
"""
function Distributions.cdf(d::InterpolatedLinear, x::Real)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]

    x <= grid[1] && return zero(x)
    x >= grid[end] && return one(x)

    return integral(itr, grid[1], x) / d.integral
end
