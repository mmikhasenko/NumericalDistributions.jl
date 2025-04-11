# NumericalDistributions.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/mmikhasenko/NumericalDistributions.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/NumericalDistributions.jl/actions)
[![Lint workflow Status](https://github.com/mmikhasenko/NumericalDistributions.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/NumericalDistributions.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mmikhasenko/NumericalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/NumericalDistributions.jl)


A Julia package for working with user-defined continuous univariate distributions where the PDF is specified numerically. The package handles normalization automatically via integration and implements sampling through numerical CDF inversion.

## Features

- Support for arbitrary continuous univariate distributions
- Automatic PDF normalization via numerical integration
- Efficient sampling through binned CDF inversion
- Support for both finite and infinite support ranges
- Implements the Distributions.jl interface

## Installation

```julia
using Pkg
Pkg.add("NumericalDistributions")
```

## Usage

Here's a simple example of creating and using a custom distribution:

```julia
using NumericalDistributions

# Define a custom distribution
const m = 0.77
const Γ = 0.15
f(x) = 1/abs2(m^2-x^2-1im*m*Γ)
dist = NumericallyIntegrable(f, (0.3, 1.5))

# The PDF is automatically normalized
pdf(dist, 0)  # returns the probability density at x=0

# Generate random samples
samples = rand(dist, 1000)

# Compute CDF
cdf(dist, 1.0)  # returns the cumulative probability at x=1.0
```

## Implementation Details

The package uses numerical integration (via `QuadGK.jl`) to normalize the PDF and compute the CDF. For sampling, it uses a binned approximation of the CDF inversion method with a linear interpolation scheme.

For distributions with infinite support, a tangent transformation is used to map the infinite interval to (-1,1) before binning.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
