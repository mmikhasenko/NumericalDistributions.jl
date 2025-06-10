"""
    invcdf(d::ContinuousUnivariateDistribution, p)

Alias for `quantile` to maintain backward compatibility. This function exists to provide compatibility with code written for earlier versions of the package.
"""
invcdf(d::ContinuousUnivariateDistribution, p) = quantile(d, p)

"""
    quantile(d::InterpolatedConstant, u::Real)

Specialized scalar quantile implementation for constant-interpolated distributions.
Delegates to the vectorized version for efficient computation using precomputed bin edges.
"""
quantile(d::InterpolatedConstant, u::Real) = quantile(d, [u])[1]

"""
    quantile(d::InterpolatedLinear, u::Real)

Specialized scalar quantile implementation for linear-interpolated distributions.
Delegates to the vectorized version which uses analytic solution to the quadratic equation
derived from the linear PDF within each bin.
"""
quantile(d::InterpolatedLinear, u::Real) = quantile(d, [u])[1]



"""
    rand(rng::AbstractRNG, d::InterpolatedConstant, n::Int64)

Fast sampling implementation for constant-interpolated distributions using the inverse transform method.
Takes advantage of the specialized `quantile` implementation to efficiently generate samples
by transforming uniform random numbers through the distribution's inverse CDF.
"""
function Base.rand(rng::AbstractRNG, d::InterpolatedConstant, n::Int64)
    u = rand(rng, n)
    return quantile(d, u)
end

"""
    rand(rng::AbstractRNG, d::InterpolatedLinear, n::Int64)

Fast sampling implementation for linearly-interpolated distributions using the inverse transform method.
Leverages the specialized `quantile` function that analytically solves the quadratic equation
arising from the inversion of the CDF with linear interpolation.
"""
function Base.rand(rng::AbstractRNG, d::InterpolatedLinear, n::Int64)
    u = rand(rng, n)
    return quantile(d, u)
end

"""
    rand(rng::AbstractRNG, d::NumericallyIntegrable, n::Int64)

Specialized sampling for general `NumericallyIntegrable` distributions with innovative handling of infinite domains.

- For finite support: Creates a constant-interpolated approximation to efficiently sample using binned CDF inversion
- For infinite support: Applies a tangent transformation to map the infinite domain to (-1,1), preserving the proper distribution shape while enabling efficient sampling

This implementation automatically adapts to the distribution's support type and ensures
accurate sampling regardless of whether the domain is bounded or unbounded.
"""
function Base.rand(rng::AbstractRNG, d::NumericallyIntegrable, n::Int64)
    if all(isfinite.(d.support))
        bD = interpolated(
            x -> d.unnormalized_pdf(x),
            range(d.support..., d.n_sampling_bins);
            degree = Constant(),
        )
        return rand(rng, bD, n)
    else
        # For infinite support, use tangent transformation
        x(z) = tan(z * π / 2)
        z(x) = atan(x) * 2 / π
        bD = interpolated(
            z -> d.unnormalized_pdf(x(z)) / cos(z * π / 2)^2,
            range(-1, 1, d.n_sampling_bins);
            degree = Constant(),
        )
        _sample = rand(rng, bD, n)
        return x.(_sample)
    end
end


# Internal methods
# precompute and broadcasted

"""
    quantile(d::InterpolatedConstant, u::AbstractVector)

Vectorized quantile computation for constant-interpolated distributions with optimized bin handling.

Implementation details:
- Constructs half-bin edges for accurate representation of the constant interpolation scheme
- Precomputes normalized weights based on bin widths for efficient lookup
- Uses binary search (`searchsortedlast`) to locate the appropriate bin for each probability
- Maps probabilities to quantiles using analytical linear transformation within each bin
"""
function quantile(d::InterpolatedConstant, u::AbstractVector)
    grid = d.unnormalized_pdf.knots[1]
    n_grid = length(grid)
    # Construct half-bin edges
    edges = Vector{eltype(grid)}(undef, n_grid + 1)
    edges[1] = grid[1]
    for i = 2:n_grid
        edges[i] = 0.5 * (grid[i] + grid[i-1])
    end
    edges[n_grid+1] = grid[end]
    values = d.unnormalized_pdf.coefs
    bin_widths = diff(edges)
    # Calculate weights
    unnormalized_weights = values .* bin_widths
    weights = unnormalized_weights ./ d.integral
    cdf_grid = cumsum(vcat(0.0, weights))
    return _invcdf_constant_scalar.(u, Ref(edges), Ref(weights), Ref(cdf_grid))
end

"""
    quantile(d::InterpolatedLinear, u::AbstractVector)

Vectorized quantile computation for linearly-interpolated distributions using quadratic equation solutions.

Implementation details:
- Normalizes PDF values from the interpolation structure
- Computes the CDF at each grid point for efficient probability-to-bin mapping
- For each probability value, locates the bin and solves the quadratic equation
  arising from the linear interpolation of the PDF within that bin
- Uses PDF values directly rather than pre-computing weights, leveraging the linearity properties
"""
function quantile(d::InterpolatedLinear, u::AbstractVector)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]
    pdf_grid = itr.coefs ./ d.integral
    cdf_grid = cdf.(Ref(d), grid)
    return _invcdf_linear_scalar.(u, Ref(grid), Ref(pdf_grid), Ref(cdf_grid))
end

"""
    _invcdf_constant_scalar(u, grid, weights, cdf_grid)

Internal helper that implements the core inversion algorithm for constant-interpolated distributions.

Key aspects of this implementation:
- Handles edge cases for u=0 and u=1 with direct mapping to support bounds
- Uses binary search with `searchsortedlast` for efficient bin location
- Implements bin-safety checks to prevent out-of-bounds access
- Computes the precise position within the bin using linear interpolation based on
  the probability position relative to the bin's cumulative probability range
"""
function _invcdf_constant_scalar(u, grid, weights, cdf_grid)
    # Handle edge cases
    u <= zero(u) && return grid[1]
    u >= one(u) && return grid[end]

    # For 0 < u < 1
    # searchsortedlast returns the index k such that cdf_grid[k] <= u < cdf_grid[k+1]
    bin_ind = searchsortedlast(cdf_grid, u)

    # Safety check for bin index
    bin_ind = max(1, min(bin_ind, length(weights)))

    # Get bin boundaries and weight
    x0, x1 = grid[bin_ind], grid[bin_ind+1]
    w = weights[bin_ind]

    # Compute position within bin
    x = x0 + (u - cdf_grid[bin_ind]) / w * (x1 - x0)
    return x
end

"""
    _invcdf_linear_scalar(u, grid, pdf_grid, cdf_grid)

Internal helper implementing quantile computation for linearly-interpolated distributions through quadratic equation solving.

Key implementation features:
- Handles boundary conditions for u=0, u=1, and edge bins
- Uses `searchsortedfirst` to locate the bin containing the target probability
- Derives and solves the quadratic equation from the linear PDF interpolation within the bin
- Selects the appropriate root based on proximity to the bin center for numerical stability
- Includes safeguards against negative discriminants by taking the maximum with zero
"""
function _invcdf_linear_scalar(u, grid, pdf_grid, cdf_grid)
    # Handle edge cases
    u <= zero(u) && return grid[1]
    u >= one(u) && return grid[end]

    # For 0 < u < 1
    idx = searchsortedfirst(cdf_grid, u)
    if idx == 1
        return grid[1]
    elseif idx > length(grid)
        return grid[end]
    else
        x0, x1 = grid[idx-1], grid[idx]
        v0, v1 = pdf_grid[idx-1], pdf_grid[idx]
        h = x1 - x0
        c0 = cdf_grid[idx-1]
        a = (v1 - v0) / h
        b = v0
        c = c0 - u

        # Solve quadratic equation for position within bin
        disc = b^2 - 2 * a * c
        disc = max(disc, 0.0)  # Ensure non-negative discriminant

        s1 = (-b + sqrt(disc)) / a
        s2 = (-b - sqrt(disc)) / a

        # Choose the solution closest to the center of the bin
        bin_center = h / 2
        s = abs(s1 - bin_center) <= abs(s2 - bin_center) ? s1 : s2

        return x0 + s
    end
end
