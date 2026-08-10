# Wrappers around GenericFFT.jl (https://github.com/JuliaApproximation/GenericFFT.jl).
# Used by the `:generic` convolution backend for element types FFTW does not support,
# including extended-precision floats. AD types (e.g. ReverseDiff.TrackedReal) are not
# supported by GenericFFT yet; `:auto` selects `:direct` for those.

using GenericFFT: generic_fft as _generic_fft_impl, generic_ifft as _generic_ifft_impl

const _generic_ifft = x -> _generic_ifft_impl(x, 1)
const _GENERIC_FFT_ELTYPE = Union{AbstractFloat, Complex{<:AbstractFloat}}

"""
    generic_fft(x::AbstractVector)

Forward FFT via [GenericFFT.jl](https://github.com/JuliaApproximation/GenericFFT.jl),
matching the `FFTW.fft` convention. Requires `eltype(x)` to be `AbstractFloat` or
`Complex{<:AbstractFloat}`.
"""
function generic_fft(x::AbstractVector{T}) where {T<:_GENERIC_FFT_ELTYPE}
    if T <: Real
        return _generic_fft_impl(complex.(x))
    end
    return _generic_fft_impl(x)
end

"""
    generic_ifft(x::AbstractVector)

Inverse FFT via GenericFFT.jl (includes `1/n` scaling).
"""
generic_ifft(x::AbstractVector{<:_GENERIC_FFT_ELTYPE}) = _generic_ifft(x)

generic_fft_supported(::Type{T}) where {T} = T <: _GENERIC_FFT_ELTYPE

generic_fft_supported(y1::AbstractVector, y2::AbstractVector) =
    generic_fft_supported(promote_type(eltype(y1), eltype(y2)))
