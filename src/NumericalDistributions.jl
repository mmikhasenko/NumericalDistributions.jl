"""
    NumericalDistributions

A Julia package for working with user-defined continuous univariate distributions
where the PDF is specified numerically. The package handles normalization automatically
via integration and implements sampling through numerical CDF inversion.

# Features
- Support for arbitrary continuous univariate distributions
- Automatic PDF normalization via numerical integration
- Efficient sampling through binned CDF inversion
- Support for both finite and infinite support ranges
"""
module NumericalDistributions

using FFTW
using Interpolations
using Distributions
using Parameters
using Random
using QuadGK

export NumericallyIntegrable

import Distributions.Statistics: mean, var
import Distributions.StatsBase: kurtosis, skewness
import Distributions: pdf, cdf, minimum, maximum, quantile

export pdf, cdf, quantile

include("types.jl")
include("moments.jl")

export integral
include("interpolate-integral.jl")

export interpolated # method
export InterpolatedLinear, InterpolatedConstant # type aliases
include("interpolated.jl")

export invcdf
include("sampling.jl")

export fft_convolve
include("convolution.jl")


end # module
