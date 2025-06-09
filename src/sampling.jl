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


function _invcdf_constant_scalar(uu, edges, values, weights, cumarr)
    uu = clamp(uu, 0.0, 1.0) # Ensure uu is strictly within [0,1]

    if uu == 0.0
        return edges[1]
    end
    if uu == 1.0
        return edges[end]
    end

    # For 0 < uu < 1
    # searchsortedlast returns the index k such that cumarr[k] <= uu < cumarr[k+1]
    # or 0 if uu < cumarr[1] (but cumarr[1] is 0, handled by uu == 0.0 check)
    binind = searchsortedlast(cumarr, uu)

    # Ensure binind is a valid index for edges and weights
    # cumarr has length ngrid+1, weights has length ngrid, edges has length ngrid+1
    # binind from searchsortedlast will be in 1:ngrid if 0 < uu < 1
    if binind < 1
        binind = 1 # Should not happen if uu > 0 and cumarr[1]=0
    elseif binind > length(weights) # length(weights) is ngrid
        binind = length(weights) # Clamp to last valid bin index for weights
    end

    σl, σr = edges[binind], edges[binind+1]
    w = weights[binind]

    if abs(w) < eps(typeof(w)) # Check if weight is effectively zero
        # CDF is flat in this region, return left edge by convention
        return σl
    end

    x = σl + (uu - cumarr[binind]) / w * (σr - σl)

    # Clamp result to the bin to handle potential floating point inaccuracies
    return clamp(x, σl, σr)
end

function fast_invcdf_constant(d::InterpolatedConstant, u::AbstractVector)
    grid = d.unnormalized_pdf.knots[1]
    ngrid = length(grid)
    # Construct half-bin edges
    edges = Vector{eltype(grid)}(undef, ngrid + 1)
    edges[1] = grid[1]
    for i = 2:ngrid
        edges[i] = 0.5 * (grid[i] + grid[i-1])
    end
    edges[ngrid+1] = grid[end]
    values = d.unnormalized_pdf.coefs
    # Bin widths: first and last are half-widths if grid is uniform
    bin_widths = diff(edges)
    # Use value at left edge for each bin (matches CDF logic)
    unnorm_weights = values .* bin_widths
    weights = unnorm_weights / d.integral
    cumarr = cumsum(vcat(0.0, weights))
    return _invcdf_constant_scalar.(u, Ref(edges), Ref(values), Ref(weights), Ref(cumarr))
end

# --- Shared fast invcdf for InterpolatedLinear ---
function _invcdf_linear_scalar(uu, grid, values, cdf_grid, integral)
    uu = clamp(uu, 0, 1)
    if uu == 0
        return grid[1]
    elseif uu == 1
        return grid[end]
    end
    idx = searchsortedfirst(cdf_grid, uu)
    if idx == 1
        return grid[1]
    elseif idx > length(grid)
        return grid[end]
    else
        x0, x1 = grid[idx-1], grid[idx]
        v0, v1 = values[idx-1], values[idx]
        h = x1 - x0
        c0 = cdf_grid[idx-1]
        a = (v1 - v0) / h / integral
        b = v0 / integral
        c = c0 - uu
        if abs(a) < 1e-14
            s = -c / b
        else
            disc = b^2 - 2 * a * c

            # Discriminant should be non-negative by the mathematical properties of CDFs
            # If it's negative, it's likely due to numerical issues - use max to avoid this
            disc = max(disc, 0.0)

            s1 = (-b + sqrt(disc)) / a
            s2 = (-b - sqrt(disc)) / a

            # Calculate distance from center of bin for both solutions
            bin_center = h / 2
            dist1 = abs(s1 - bin_center)
            dist2 = abs(s2 - bin_center)

            # Select the solution closest to the center of the bin
            # and clamp to valid range if needed
            s = dist1 <= dist2 ? s1 : s2

            # Clamp to ensure we're in the valid range [0,h]
            s = clamp(s, 0, h)
        end
        return x0 + s
    end
end

function fast_invcdf_linear(d::InterpolatedLinear, u::AbstractVector)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]
    values = itr.coefs
    cdf_grid = cdf.(Ref(d), grid)
    integral = d.integral
    return _invcdf_linear_scalar.(u, Ref(grid), Ref(values), Ref(cdf_grid), Ref(integral))
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
