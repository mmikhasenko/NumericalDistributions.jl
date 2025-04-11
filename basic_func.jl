begin
    using Distributions
    using Random
    using QuadGK
end

begin
    struct binned1dDensity{T}
        grid::Vector{Float64}
        cumarr::Vector{T}
    end
    function getbinned1dDensity(g, lims, nBins)
        grid = collect(range(lims..., nBins))
        bin_centers = (grid[2:end] .+ grid[1:end-1]) ./ 2
        unnorm_weights = map(g, bin_centers)
        weights = unnorm_weights ./ sum(unnorm_weights)
        cumarr = pushfirst!(cumsum(vcat(weights...), dims = 1), 0)
        return binned1dDensity(grid, cumarr)
    end
    function Base.rand(rng::AbstractRNG, bD::binned1dDensity)
        binind = findfirst(bD.cumarr .> rand()) - 1
        σl, σr = bD.grid[binind], bD.grid[binind+1]
        σ = σl + rand(rng) * (σr - σl)
        return σ
    end
end

begin
    struct NumericallyIntegrable{F, S} <: ContinuousUnivariateDistribution
        unnormalized_pdf::F
        integral::Float64
        support::S
        n_sampling_bins::Int
        function NumericallyIntegrable(f, support = (-Inf, Inf), n_sampling_bins = 300)
            integral = quadgk(f, support...)[1]
            F, S = typeof(f), typeof(support)
            new{F, S}(f, integral, support, n_sampling_bins)
        end
    end
    function Distributions.pdf(d::NumericallyIntegrable, m::Real)
        return d.support[1] < m < d.support[2] ?
               d.unnormalized_pdf(m) / d.integral :
               zero(m)
    end
    function Distributions.cdf(d::NumericallyIntegrable, m::Real)
        m < support[1] && return zero(m)
        m > support[2] && return one(m)
        quadgk(d.unnormalized_pdf, support[1], m)[1] / d.integral
    end
    function Base.rand(rng::AbstractRNG, d::NumericallyIntegrable, n::Int64 = 1)
        if all(isfinite.(d.support))
            d = getbinned1dDensity(d.unnormalized_pdf, d.support, d.n_sampling_bins)
            return map(1:n) do _
                rand(rng, d)
            end
        end
        # infinite support range - use mapping
        x(z) = tan(z * π / 2)
        z(x) = atan(x) / π * 2
        d = getbinned1dDensity((-1, 1), d.n_sampling_bins) do z
            d.unnormalized_pdf(x(z)) / cos(z * π / 2)^2
        end
        _sample = map(1:n) do _
            rand(rng, d)
        end
        return x.(_sample)
    end
end
