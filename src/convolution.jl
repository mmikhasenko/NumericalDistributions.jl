"""
    _fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0, pow2=true)

Internal: Numerically convolve two vectors of values using FFTW. Returns the new grid and convolution result as vectors.
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
    fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0, pow2=true, degree=Linear())

Numerically convolve two vectors of values and return a NumericallyIntegrable distribution.
"""
function fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1, t0_2, pow2 = true)
    t, h = _fft_convolve(y1, y2; Δ, t0_1, t0_2, pow2)
    itr = interpolate((t,), h, Gridded(Linear()))
    return NumericallyIntegrable(itr, extrema(t); n_sampling_bins = length(t))
end

"""
    fft_convolve(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution; pow2=true, degree=Linear(), gridsize=1000)

Numerically convolve two probability density functions (PDFs) using FFTW. Accepts any objects that are subtypes of `ContinuousUnivariateDistribution` (e.g., `NumericallyIntegrable`, `Interpolated`, or standard continuous distributions).
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
