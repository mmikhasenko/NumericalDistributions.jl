"""
    _fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0, pow2=true)

Internal FFT-based convolution implementation optimized for performance.
Key optimizations include:
- Optional power-of-2 padding for faster FFT operations
- In-place zero padding to minimize memory allocation
- Domain tracking to correctly position the resulting convolution
- Direct use of FFTW.jl for efficient computation
"""
function _fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1, t0_2, pow2 = true)
    M, N = length(y1), length(y2)
    Lraw = M + N - 1
    L = pow2 ? nextpow(2, Lraw) : Lraw
    F = fft([y1; zeros(L - M)])
    G = fft([y2; zeros(L - N)])
    h_full = real(ifft(F .* G)) * Δ
    h = h_full[1:Lraw]
    t0_h = t0_1 + t0_2
    t = t0_h .+ (0:Lraw-1) .* Δ
    return t, h
end

"""
    fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0, pow2=true)

FFT-accelerated convolution for raw PDF vectors that preserves normalization and domain information.
This lower-level interface directly accepts PDF values on a grid and handles the grid properly
by tracking the domain boundaries. The result is automatically wrapped as a distribution object
with linear interpolation for smooth representation.
"""
function fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1, t0_2, pow2 = true)
    t, h = _fft_convolve(y1, y2; Δ, t0_1, t0_2, pow2)
    itr = interpolate((t,), h, Gridded(Linear()))
    return NumericallyIntegrable(itr, extrema(t); n_sampling_bins = length(t))
end

"""
    fft_convolve(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution;
                pow2=true, gridsize=1000)

Computes the distribution of the convolution of two independent random variables using FFT.

The method works with univariate distributions (`<:Distributions.ContinuousUnivariateDistribution`).
The grid spacing is determined using the supports of the input distributions and gridsize.

# Arguments
- `d1::ContinuousUnivariateDistribution`: First distribution to convolve
- `d2::ContinuousUnivariateDistribution`: Second distribution to convolve
- `pow2::Bool=true`: Whether to pad arrays to powers of 2 for faster FFT computation
- `gridsize::Int=1000`: Number of grid points to use when sampling each distribution.
  Higher values increase accuracy but also increase memory usage and computation time

# Returns
- `NumericallyIntegrable`: A distribution representing the convolution of `d1` and `d2`,
  implemented as a linearly interpolated numerical distribution

# Notes
- Both input distributions must have finite support. For distributions with infinite support,
  use `Truncated` from Distributions.jl to create finite support approximations
- The resulting distribution's support will be the sum of the supports of the input distributions
- Computational complexity is O(n log n) where n is proportional to gridsize

# Example
````julia
using NumericalDistributions
using Distributions

d1 = truncated(Normal(0, 1), -3, 3)
d2 = truncated(Gamma(2, 1), 0, 10)
d_conv = fft_convolve(d1, d2)
````

"""
function fft_convolve(
    d1::ContinuousUnivariateDistribution,
    d2::ContinuousUnivariateDistribution;
    pow2 = true,
    gridsize = 1000,
)
    a1, b1 = minimum(d1), maximum(d1)
    a2, b2 = minimum(d2), maximum(d2)
    if !isfinite(a1) || !isfinite(b1) || !isfinite(a2) || !isfinite(b2)
        error(
            "fft_convolve: Both distributions must have finite support. Consider truncating the support using a suitable interval, e.g., `Truncated(d, -5, 5)` for a Normal distribution.",
        )
    end
    Δ1 = (b1 - a1) / (gridsize - 1)
    Δ2 = (b2 - a2) / (gridsize - 1)
    Δ = min(Δ1, Δ2)
    total_steps = length(a1:Δ:b1) + length(a2:Δ:b2) - 1
    if total_steps > 2^16
        @warn """
The total number of grid points ($(total_steps)) is very large.
It is computed as `min(Δ1, Δ2)`, where `Δ1` and `Δ2` are the step sizes for the two distributions.
    Δ1 = $(Δ1)
    Δ2 = $(Δ2)
Consider reducing `gridsize`, or restricting the support of your distributions.
"""
    end
    t1 = a1:Δ:b1
    t2 = a2:Δ:b2
    y1 = pdf.(d1, t1)
    y2 = pdf.(d2, t2)
    return fft_convolve(y1, y2; Δ, t0_1 = a1, t0_2 = a2, pow2)
end
