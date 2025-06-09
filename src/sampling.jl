"""
    invcdf(d::InterpolatedConstant, u::Real)

Efficiently invert the CDF for an InterpolatedConstant distribution.
Finds the bin for u and computes the analytic location within the bin, accounting for the value change in the middle of the bin.
"""
invcdf(d::InterpolatedConstant, u::Real) = fast_invcdf_constant(d, [u])[1]

"""
    invcdf(d::InterpolatedLinear, u::Real)

Efficiently invert the CDF for an InterpolatedLinear distribution.
Finds the bin for u and computes the analytic location within the bin using the linear CDF formula.
"""
invcdf(d::InterpolatedLinear, u::Real) = fast_invcdf_linear(d, [u])[1]


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

function fast_invcdf_constant(d::InterpolatedConstant, u::AbstractVector)
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

function fast_invcdf_linear(d::InterpolatedLinear, u::AbstractVector)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]
    pdf_grid = itr.coefs ./ d.integral
    cdf_grid = cdf.(Ref(d), grid)
    return _invcdf_linear_scalar.(u, Ref(grid), Ref(pdf_grid), Ref(cdf_grid))
end

# Fast sampling for InterpolatedConstant using shared invcdf
function Base.rand(rng::AbstractRNG, d::InterpolatedConstant, n::Int64 = 1)
    u = rand(rng, n)
    return fast_invcdf_constant(d, u)
end

# Fast sampling for InterpolatedLinear using shared invcdf
function Base.rand(rng::AbstractRNG, d::InterpolatedLinear, n::Int64 = 1)
    u = rand(rng, n)
    return fast_invcdf_linear(d, u)
end

Base.rand(d::NumericallyIntegrable, n::Int64 = 1) = rand(Random.GLOBAL_RNG, d, n)
