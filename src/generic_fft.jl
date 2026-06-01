# Parametric FFT/IFFT (Cooley–Tukey for power-of-two lengths, naive DFT otherwise).
# Intended for types that FFTW does not support (TrackedReal, Dual, etc.).
# A future GenericFFT.jl package can supersede this module.

"""
    generic_fft(x::AbstractVector)

Discrete Fourier transform using a type-parametric algorithm (`+`, `*`, `cis`).
Matches the `FFTW.fft` convention (forward sign ``-2\\pi i k / n``).

For length `n` that is a power of two, uses Cooley–Tukey ``O(n \\log n)``.
Otherwise uses a naive ``O(n^2)`` DFT (suitable for small `n` only).
"""
function generic_fft(x::AbstractVector)
    n = length(x)
    n == 0 && return eltype(x)[]
    Tx = eltype(x)
    Tc = typeof(complex(one(Tx)))
    xc = map(Tc, x)
    if ispow2(n)
        return _generic_fft_pow2(xc)
    end
    return _generic_dft_forward(xc)
end

"""
    generic_ifft(x::AbstractVector)

Inverse DFT matching `FFTW.ifft` (includes `1/n` scaling).
"""
function generic_ifft(x::AbstractVector)
    n = length(x)
    n == 0 && return eltype(x)[]
    y = conj.(generic_fft(conj.(x)))
    invn = one(eltype(y)) / n
    return y .* invn
end

function _generic_fft_pow2(x::AbstractVector{T}) where {T}
    n = length(x)
    n == 1 && return copy(x)
    @assert ispow2(n)
    half = n ÷ 2
    ev = _generic_fft_pow2(@view x[1:2:n-1])
    od = _generic_fft_pow2(@view x[2:2:n])
    out = similar(x, n)
    ω = cis(-2 * pi / n)
    w = one(T)
    @inbounds for k in 1:half
        t = w * od[k]
        ek = ev[k]
        out[k] = ek + t
        out[k + half] = ek - t
        w *= ω
    end
    return out
end

function _generic_dft_forward(x::AbstractVector{T}) where {T}
    n = length(x)
    out = similar(x, n)
    @inbounds for k in 1:n
        s = zero(T)
        for j in 1:n
            s += x[j] * cis(-2 * pi * (k - 1) * (j - 1) / n)
        end
        out[k] = s
    end
    return out
end
