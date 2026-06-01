# Convolution backends for sampled PDFs on a uniform grid.
#
# Backends (`algorithm` keyword):
#   :fftw    — FFTW via AbstractFFTs (Float32/64 and complex counterparts)
#   :generic — parametric FFT in generic_fft.jl (AD-friendly types)
#   :direct  — O(M·N) reference convolution
#   :auto    — :fftw when supported, else :generic (or :direct if padded length is huge)

const _FFTW_CONVOLVE_ELTYPE = Union{Float32, Float64, Complex{Float32}, Complex{Float64}}

const _CONVOLUTION_ALGORITHMS = (:auto, :fftw, :generic, :direct)
# Deprecated alias kept for compatibility with earlier PR draft
const _CONVOLUTION_ALGORITHM_ALIASES = Dict(:fft => :fftw)

_fftw_convolve_supported(y1::AbstractVector, y2::AbstractVector) =
    promote_type(eltype(y1), eltype(y2)) <: _FFTW_CONVOLVE_ELTYPE

"""Maximum padded FFT length for `:auto` to select `:generic` (else `:direct`)."""
const GENERIC_FFT_AUTO_MAX_LENGTH = 16_384

function _normalize_convolution_algorithm(algorithm::Symbol)
    algorithm = get(_CONVOLUTION_ALGORITHM_ALIASES, algorithm, algorithm)
    algorithm in _CONVOLUTION_ALGORITHMS ||
        error(
            "fft_convolve: `algorithm` must be one of $(collect(_CONVOLUTION_ALGORITHMS)) (got $(repr(algorithm)))",
        )
    return algorithm
end

function _padded_convolution_length(y1::AbstractVector, y2::AbstractVector; pow2::Bool)
    Lraw = length(y1) + length(y2) - 1
    return pow2 ? nextpow(2, Lraw) : Lraw, Lraw
end

function _default_convolution_algorithm(y1::AbstractVector, y2::AbstractVector; pow2::Bool)
    if _fftw_convolve_supported(y1, y2)
        return :fftw
    end
    L, _ = _padded_convolution_length(y1, y2; pow2)
    return L <= GENERIC_FFT_AUTO_MAX_LENGTH ? :generic : :direct
end

function _resolve_convolution_algorithm(
    y1::AbstractVector,
    y2::AbstractVector;
    algorithm::Symbol,
    pow2::Bool,
)
    algorithm = _normalize_convolution_algorithm(algorithm)
    if algorithm === :auto
        return _default_convolution_algorithm(y1, y2; pow2)
    end
    if algorithm === :fftw && !_fftw_convolve_supported(y1, y2)
        T = promote_type(eltype(y1), eltype(y2))
        error(
            "fft_convolve: FFTW backend does not support element type $T; use `algorithm = :generic`, `:direct`, or `:auto`",
        )
    end
    return algorithm
end

"""
    _direct_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0)

O(M·N) discrete convolution. AD-safe reference implementation.
"""
function _direct_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1, t0_2)
    M, N = length(y1), length(y2)
    Lraw = M + N - 1
    T = promote_type(eltype(y1), eltype(y2), typeof(Δ))
    h = zeros(T, Lraw)
    @inbounds for i in 1:M
        yi = y1[i]
        for j in 1:N
            h[i + j - 1] += yi * y2[j]
        end
    end
    h .*= Δ
    t0_h = t0_1 + t0_2
    t = t0_h .+ (0:Lraw - 1) .* Δ
    return t, h
end

function _fft_convolve_impl(
    y1::AbstractVector,
    y2::AbstractVector,
    fft_fn,
    ifft_fn;
    Δ,
    t0_1,
    t0_2,
    pow2 = true,
)
    M, N = length(y1), length(y2)
    Lraw = M + N - 1
    L = pow2 ? nextpow(2, Lraw) : Lraw
    T = promote_type(eltype(y1), eltype(y2), typeof(Δ))
    z = zero(T)
    F = fft_fn([y1; fill(z, L - M)])
    G = fft_fn([y2; fill(z, L - N)])
    h_full = real.(ifft_fn(F .* G)) .* Δ
    h = h_full[1:Lraw]
    t0_h = t0_1 + t0_2
    t = t0_h .+ (0:Lraw - 1) .* Δ
    return t, h
end

"""FFTW-backed convolution (fast path for machine floats)."""
_fft_convolve(y1::AbstractVector, y2::AbstractVector; kwargs...) =
    _fft_convolve_impl(y1, y2, fft, ifft; kwargs...)

"""Parametric-FFT convolution (traceable element types)."""
_generic_fft_convolve(y1::AbstractVector, y2::AbstractVector; kwargs...) =
    _fft_convolve_impl(y1, y2, generic_fft, generic_ifft; kwargs...)

function _convolve_vectors(
    y1::AbstractVector,
    y2::AbstractVector;
    Δ,
    t0_1,
    t0_2,
    pow2 = true,
    algorithm = :auto,
)
    alg = _resolve_convolution_algorithm(y1, y2; algorithm, pow2)
    if alg === :fftw
        return _fft_convolve(y1, y2; Δ, t0_1, t0_2, pow2)
    elseif alg === :generic
        return _generic_fft_convolve(y1, y2; Δ, t0_1, t0_2, pow2)
    elseif alg === :direct
        return _direct_convolve(y1, y2; Δ, t0_1, t0_2)
    else
        error("fft_convolve: internal error, unresolved algorithm $alg")
    end
end

"""
    fft_convolve(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0,
                 pow2=true, algorithm=:auto)

Convolve two sampled PDF vectors on uniform grids with spacing `Δ`, returning a
`NumericallyIntegrable` distribution on the convolution support.

# Convolution backends (`algorithm`)
- `:auto` — `:fftw` for `Float32`/`Float64`/complex machine types; otherwise `:generic`
  when the padded length is ≤ `GENERIC_FFT_AUTO_MAX_LENGTH`, else `:direct`
- `:fftw` — FFTW (requires FFTW-supported element types)
- `:generic` — type-parametric FFT (e.g. `ReverseDiff.TrackedReal`, `ForwardDiff.Dual`)
- `:direct` — O(M·N) reference convolution (always AD-traceable, slow for large grids)

`pow2` pads to the next power of two before FFT-based backends (ignored by `:direct`).

See also [`generic_fft`](@ref), [`generic_ifft`](@ref).
"""
function fft_convolve(
    y1::AbstractVector,
    y2::AbstractVector;
    Δ,
    t0_1,
    t0_2,
    pow2 = true,
    algorithm = :auto,
)
    t, h = _convolve_vectors(y1, y2; Δ, t0_1, t0_2, pow2, algorithm)
    itr = interpolate((t,), h, Gridded(Linear()))
    return NumericallyIntegrable(itr, extrema(t); n_sampling_bins = length(t))
end

"""
    fft_convolve(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution;
                pow2=true, gridsize=1000, algorithm=:auto)

Computes the distribution of the convolution of two independent random variables.

See the vector method [`fft_convolve`](@ref) for the `algorithm` keyword.
"""
function fft_convolve(
    d1::ContinuousUnivariateDistribution,
    d2::ContinuousUnivariateDistribution;
    pow2 = true,
    gridsize = 1000,
    algorithm = :auto,
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
    return fft_convolve(y1, y2; Δ, t0_1 = a1, t0_2 = a2, pow2, algorithm)
end
