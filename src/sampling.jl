"""
    invcdf(d::InterpolatedConstant, u::Real)

Efficiently invert the CDF for an InterpolatedConstant distribution.
Finds the bin for u and computes the analytic location within the bin, accounting for the value change in the middle of the bin.
"""
function invcdf(d::InterpolatedConstant, u::Real)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]
    ngrid = length(grid)
    # Reconstruct bin edges as in cdf
    edges = Vector{eltype(grid)}(undef, ngrid + 1)
    edges[1] = grid[1]
    for i = 2:ngrid
        edges[i] = 0.5 * (grid[i] + grid[i-1])
    end
    edges[ngrid+1] = grid[end]
    cdf_edges = cdf.(Ref(d), edges)
    # Clamp u to [0,1]
    u = clamp(u, 0, 1)
    idx = searchsortedfirst(cdf_edges, u)
    if idx == 1
        return edges[1]
    elseif idx > length(edges)
        return edges[end]
    else
        # Uniform sampling within the bin
        x0, x1 = edges[idx-1], edges[idx]
        c0, c1 = cdf_edges[idx-1], cdf_edges[idx]
        t = (u - c0) / (c1 - c0)
        return x0 + t * (x1 - x0)
    end
end

"""
    invcdf(d::InterpolatedLinear, u::Real)

Efficiently invert the CDF for an InterpolatedLinear distribution.
Finds the bin for u and computes the analytic location within the bin using the linear CDF formula.
"""
function invcdf(d::InterpolatedLinear, u::Real)
    itr = d.unnormalized_pdf
    grid = itr.knots[1]
    # Handle edge cases
    u == zero(u) && return grid[1]
    u == one(u) && return grid[end]
    # Handle general case
    values = itr.coefs
    n = length(grid)
    cdf_grid = cdf.(Ref(d), grid)
    u = clamp(u, 0, 1)
    idx = searchsortedfirst(cdf_grid, u)
    if idx == 1
        return grid[1]
    elseif idx > length(grid)
        return grid[end]
    else
        x0, x1 = grid[idx-1], grid[idx]
        v0, v1 = values[idx-1], values[idx]
        h = x1 - x0
        c0 = cdf_grid[idx-1]
        a = (v1 - v0) / h / d.integral
        b = v0 / d.integral
        c = c0 - u
        if abs(a) < 1e-14
            # Degenerate to linear
            s = -c / b
        else
            disc = b^2 - 2 * a * c
            @assert disc >= 0 "Negative discriminant in invcdf for InterpolatedLinear."
            s1 = (-b + sqrt(disc)) / a
            s2 = (-b - sqrt(disc)) / a
            # Select s in [0, h]
            s_candidates = filter(s -> 0 <= s <= h, (s1, s2))
            @assert !isempty(s_candidates) "No valid solution for invcdf in bin."
            s = s_candidates[1]
        end
        return x0 + s
    end
end


Base.rand(rng::AbstractRNG, d::InterpolatedConstant, n::Int64 = 1) =
    invcdf.(Ref(d), rand(rng, n))

Base.rand(rng::AbstractRNG, d::InterpolatedLinear, n::Int64 = 1) =
    invcdf.(Ref(d), rand(rng, n))

function Base.rand(rng::AbstractRNG, d::NumericallyIntegrable, n::Int64 = 1)
    if all(isfinite.(d.support))
        bD = BinnedDensity(x -> d.unnormalized_pdf(x), d.support, d.n_sampling_bins)
        return map(_ -> rand(rng, bD), 1:n)
    else
        # For infinite support, use tangent transformation
        x(z) = tan(z * π / 2)
        z(x) = atan(x) * 2 / π
        bD = BinnedDensity(
            z -> d.unnormalized_pdf(x(z)) / cos(z * π / 2)^2,
            (-1, 1),
            d.n_sampling_bins,
        )
        _sample = map(_ -> rand(rng, bD), 1:n)
        return x.(_sample)
    end
end



"""
    BinnedDensity{T}

A type representing a binned approximation of a 1D probability density function.
Used internally for efficient sampling from NumericallyIntegrable distributions.

# Fields
- `grid::Vector{Float64}`: The bin edges
- `cumarr::Vector{T}`: The cumulative probabilities at bin edges
"""
struct BinnedDensity{T}
    grid::Vector{Float64}
    cumarr::Vector{T}
end

"""
    BinnedDensity(g, lims, nBins)

Create a binned representation of a 1D probability density function.

# Arguments
- `g`: The probability density function
- `lims`: Tuple of (lower, upper) limits for binning
- `nBins`: Number of bins to use

# Returns
- `BinnedDensity`: A binned representation of the PDF
"""
function BinnedDensity(g, lims, nBins)
    # Create nBins + 1 points to have the correct number of bin edges
    grid = collect(range(lims..., length = nBins + 1))
    bin_centers = (grid[2:end] .+ grid[1:end-1]) ./ 2
    unnorm_weights = map(g, bin_centers)
    weights = unnorm_weights ./ sum(unnorm_weights)
    cumarr = pushfirst!(cumsum(vcat(weights...), dims = 1), 0)
    return BinnedDensity(grid, cumarr)
end

function Base.rand(rng::AbstractRNG, bD::BinnedDensity)
    binind = findfirst(bD.cumarr .> rand(rng)) - 1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl + rand(rng) * (σr - σl)
    return σ
end


Base.rand(d::NumericallyIntegrable, n::Int64 = 1) = rand(Random.GLOBAL_RNG, d, n)
