"""
    binned1dDensity{T}

A type representing a binned approximation of a 1D probability density function.
Used internally for efficient sampling from NumericallyIntegrable distributions.

# Fields
- `grid::Vector{Float64}`: The bin edges
- `cumarr::Vector{T}`: The cumulative probabilities at bin edges
"""
struct binned1dDensity{T}
    grid::Vector{Float64}
    cumarr::Vector{T}
end

"""
    getbinned1dDensity(g, lims, nBins)

Create a binned representation of a 1D probability density function.

# Arguments
- `g`: The probability density function
- `lims`: Tuple of (lower, upper) limits for binning
- `nBins`: Number of bins to use

# Returns
- `binned1dDensity`: A binned representation of the PDF
"""
function getbinned1dDensity(g, lims, nBins)
    grid = collect(range(lims..., nBins))
    bin_centers = (grid[2:end] .+ grid[1:end-1]) ./ 2
    unnorm_weights = map(g, bin_centers)
    weights = unnorm_weights ./ sum(unnorm_weights)
    cumarr = pushfirst!(cumsum(vcat(weights...), dims = 1), 0)
    return binned1dDensity(grid, cumarr)
end

function Base.rand(rng::AbstractRNG, bD::binned1dDensity)
    binind = findfirst(bD.cumarr .> rand(rng)) - 1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl + rand(rng) * (σr - σl)
    return σ
end

"""
    rand(rng::AbstractRNG, d::NumericallyIntegrable[, n::Int64=1])

Generate random samples from a NumericallyIntegrable distribution.

For finite support, uses direct binning of the PDF. For infinite support,
applies a coordinate transformation to map the infinite interval to (-1,1)
before binning.
"""
function Base.rand(rng::AbstractRNG, d::NumericallyIntegrable, n::Int64 = 1)
    if all(isfinite.(d.support))
        bD = getbinned1dDensity(
            x -> d.unnormalized_pdf(x),
            d.support,
            d.n_sampling_bins)
        return map(_ -> rand(rng, bD), 1:n)
    end

    # For infinite support, use tangent transformation
    x(z) = tan(z * π / 2)
    z(x) = atan(x) * 2 / π
    bD = getbinned1dDensity(
        z -> d.unnormalized_pdf(x(z)) / cos(z * π / 2)^2,
        (-1, 1),
        d.n_sampling_bins)
    _sample = map(_ -> rand(rng, bD), 1:n)
    return x.(_sample)
end

Base.rand(d::NumericallyIntegrable, n::Int64 = 1) = rand(Random.GLOBAL_RNG, d, n)
