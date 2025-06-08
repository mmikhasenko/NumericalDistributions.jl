"""
    Interpolated(f, binedges; degree = Linear())

Create a distribution using interpolation of a function over a grid of points.

# Arguments
- `f`: The function to be interpolated
- `binedges::AbstractVector`: Grid points for interpolation
- `degree`: Interpolation type, either `Linear()` or `Constant()`

# Returns
- `NumericallyIntegrable`: A distribution based on the interpolated function
"""
function Interpolated(f, binedges::AbstractVector; degree = Linear())
    yv_unnorm = f.(binedges)
    NumericallyIntegrable(
        interpolate((binedges,), yv_unnorm, Gridded(degree)),
        extrema(binedges),
        n_sampling_bins = length(binedges),
    )
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
    values = itr.coefs
    n = length(grid)

    # Construct bin edges - these are not equally spaced
    edges = Vector{eltype(grid)}(undef, n + 1)
    edges[1] = grid[1]
    for i = 2:n
        edges[i] = 0.5 * (grid[i] + grid[i-1])
    end
    edges[n+1] = grid[n]

    # Handle out-of-bounds
    x <= edges[1] && return zero(x)
    x >= edges[end] && return one(x)

    idx = searchsortedlast(edges, x)
    s = zero(x)
    for i = 1:idx-1
        s += (edges[i+1] - edges[i]) * values[i]
    end
    s += (x - edges[idx]) * values[idx]

    return s / d.integral
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
    values = itr.coefs
    n = length(grid)

    x <= grid[1] && return zero(x)
    x >= grid[end] && return one(x)

    idx = searchsortedlast(grid, x)
    s = zero(x)

    # Sum over full bins before x
    for i = 1:idx-1
        h = grid[i+1] - grid[i]
        v0, v1 = values[i], values[i+1]
        s += h * v0 + 0.5 * (v1 - v0) * h
    end

    # Partial bin
    δ = x - grid[idx]
    h = grid[idx+1] - grid[idx]
    v0, v1 = values[idx], values[idx+1]
    s += δ * v0 + 0.5 * (v1 - v0) / h * δ^2

    return s / d.integral
end
