

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

