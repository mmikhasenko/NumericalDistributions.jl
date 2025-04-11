# NumericalDistributions.jl

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

# Define a custom distribution (truncated normal distribution)
f(x) = exp(-x^2/2) * (abs(x) < 2)
dist = NumericallyIntegrable(f, (-2, 2))

# The PDF is automatically normalized
pdf(dist, 0)  # returns the probability density at x=0

# Generate random samples
samples = rand(dist, 1000)

# Compute CDF
cdf(dist, 1.0)  # returns the cumulative probability at x=1.0
```

## Implementation Details

The package uses numerical integration (via QuadGK.jl) to normalize the PDF and compute the CDF. For sampling, it uses a binned approximation of the CDF inversion method:

1. The support range is divided into bins
2. The PDF is evaluated at bin centers and normalized
3. The CDF is approximated by cumulative sum of bin probabilities
4. Random samples are generated by:
   - Drawing a uniform random number u ∈ [0,1]
   - Finding the bin containing u in the cumulative probabilities
   - Linear interpolation within the bin

For distributions with infinite support, a tangent transformation is used to map the infinite interval to (-1,1) before binning.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.