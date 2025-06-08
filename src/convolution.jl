"""
    convolve_vectors(y1::AbstractVector, y2::AbstractVector; Δ, t0_1=0.0, t0_2=0.0, pow2=true, degree=Linear())

Numerically convolve two vectors of values (not necessarily PDFs) defined on the same grid using FFTW.

# Arguments
- `y1`, `y2`: Vectors of values (e.g., sampled PDF values).
- `Δ`: Step size of the grid (required).
- `t0_1`, `t0_2`: Starting points of the grids.
- `pow2`: If true, zero-pad to the next power of 2 for FFT efficiency.
- `degree`: Interpolation type for the result (default: `Linear()`).

# Returns
- An `Interpolated` object representing the convolution result on the extended grid.
"""
function convolve_vectors(
    y1::AbstractVector,
    y2::AbstractVector;
    Δ,
    t0_1 = 0.0,
    t0_2 = 0.0,
    pow2 = true,
    degree = Linear(),
)
    M, N = length(y1), length(y2)
    Lraw = M + N - 1
    L = pow2 ? nextpow(2, Lraw) : Lraw
    F = fft([y1; zeros(L - M)])
    G = fft([y2; zeros(L - N)])
    h_full = real(ifft(F .* G)) * Δ
    h = h_full[1:Lraw]
    t0_h = t0_1 + t0_2
    t = t0_h .+ (0:Lraw-1) .* Δ
    itr = interpolate((t,), h, Gridded(degree))
    return NumericallyIntegrable(itr, extrema(t), length(t))
end

"""
    convolve_pdfs(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution; pow2=true, degree=Linear(), gridsize=1000)

Numerically convolve two probability density functions (PDFs) using FFTW. Accepts any objects that are subtypes of `ContinuousUnivariateDistribution` (e.g., `NumericallyIntegrable`, `Interpolated`, or standard continuous distributions).

# Arguments
- `d1`, `d2`: Distributions (subtypes of `ContinuousUnivariateDistribution`).
- `pow2`: If true, zero-pad to the next power of 2 for FFT efficiency.
- `degree`: Interpolation type for the result (default: `Linear()`).
- `gridsize`: Number of points to sample each PDF on the grid.

# Returns
- An `Interpolated` object representing the convolution result on the extended grid.

# Warning
Both `d1` and `d2` must be normalized probability density functions (i.e., their integrals must be 1).
"""
function convolve_pdfs(
    d1::ContinuousUnivariateDistribution,
    d2::ContinuousUnivariateDistribution;
    pow2 = true,
    degree = Linear(),
    gridsize = 1000,
)
    a1, b1 = minimum(d1), maximum(d1)
    a2, b2 = minimum(d2), maximum(d2)
    # Check for infinite support and warn or restrict
    if !isfinite(a1) || !isfinite(b1) || !isfinite(a2) || !isfinite(b2)
        error(
            "convolve_pdfs: Both distributions must have finite support. Consider truncating the support using a suitable interval, e.g., `Truncated(d, -5, 5)` for a Normal distribution, `d = Normal()`.",
        )
    end
    Δ1 = (b1 - a1) / (gridsize - 1)
    Δ2 = (b2 - a2) / (gridsize - 1)
    Δ = min(Δ1, Δ2)
    # Warn if the number of steps is very large
    total_steps = length(a1:Δ:b1) + length(a2:Δ:b2) - 1
    if total_steps > 2^16  # 65536
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
    return convolve_vectors(
        y1,
        y2;
        Δ = Δ,
        t0_1 = a1,
        t0_2 = a2,
        pow2 = pow2,
        degree = degree,
    )
end

# Update InterpolatedLinear/Constant methods to call convolve_pdfs
default_interpolated_call(d1, d2; pow2 = true, degree = Linear(), gridsize = 1000) =
    convolve_pdfs(d1, d2; pow2 = pow2, degree = degree, gridsize = gridsize)

function convolve_pdfs(
    d1::InterpolatedLinear,
    d2::InterpolatedLinear;
    pow2 = true,
    degree = Linear(),
    gridsize = 1000,
)
    default_interpolated_call(d1, d2; pow2 = pow2, degree = degree, gridsize = gridsize)
end

function convolve_pdfs(
    d1::InterpolatedConstant,
    d2::InterpolatedConstant;
    pow2 = true,
    degree = Constant(),
    gridsize = 1000,
)
    default_interpolated_call(d1, d2; pow2 = pow2, degree = degree, gridsize = gridsize)
end
