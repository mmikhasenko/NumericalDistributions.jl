"""
    NumericallyInvertibleDistributions

A Julia package for working with user-defined continuous univariate distributions
where the PDF is specified numerically. The package handles normalization automatically
via integration and implements sampling through numerical CDF inversion.

# Features
- Support for arbitrary continuous univariate distributions
- Automatic PDF normalization via numerical integration
- Efficient sampling through binned CDF inversion
- Support for both finite and infinite support ranges
"""
module NumericallyInvertibleDistributions

using Distributions
using Random
using QuadGK

export NumericallyIntegrable
export binned1dDensity, getbinned1dDensity

include("types.jl")
include("sampling.jl")

end # module
