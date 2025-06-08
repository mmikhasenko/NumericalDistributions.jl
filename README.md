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
- Interpolated distributions with constant or linear degree options
- Implements the Distributions.jl interface

## Installation

```julia
julia> ]
pkg> add https://github.com/mmikhasenko/NumericalDistributions.jl
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
pdf(dist, 0.77)  # returns the probability density at x=0.77

# Generate random samples
samples = rand(dist, 1000)

# Compute CDF
cdf(dist, 1.0)  # returns the cumulative probability at x=1.0
```

### Interpolated Distributions

The package provides an `Interpolated` constructor for creating distributions based on interpolated functions:

```julia
using NumericalDistributions

# Create an interpolated distribution with constant degree
a = Interpolated(x -> abs(x), -2:0.5:1; degree = Constant())
pdf(a, 0.4)  # Evaluates PDF at x=0.4
cdf(a, 0.0)  # Computes CDF at x=0.0

# Create an interpolated distribution with linear degree
b = Interpolated(x -> abs(x), -2:0.5:1; degree = Linear())
pdf(b, 0.4)  # Evaluates PDF with linear interpolation
cdf(b, 0.0)  # Computes CDF with linear interpolation
```

## Numerical Convolution of PDFs (via FFTW)

NumericalDistributions.jl supports fast convolution of user-defined, numerically specified probability density functions (PDFs) using the Fast Fourier Transform (FFTW). This allows you to compute the PDF of the sum of two independent random variables, each with a numerically defined distribution.

### Available Functions

- `fft_convolve(d1, d2; ...)`: Convolve any two distributions (subtypes of `ContinuousUnivariateDistribution`). Automatically samples both on a common grid and returns the convolution as a NumericallyIntegrable distribution.

### Example: Convolution of Two PDFs

```julia
using NumericalDistributions

# Define two custom PDFs on different supports
x1 = -2:0.01:2
x2 = -1:0.01:8
f1(x) = exp(-x^4)                    # Even, super-Gaussian
f2(x) = 1 / ((x - 3)^2 + 1)          # Cauchy-like, centered at 3

# Create interpolated distributions
A = Interpolated(f1, x1)
B = Interpolated(f2, x2)

# Convolve the two distributions
C = fft_convolve(A, B)

# Evaluate the resulting PDF at a point
pdf(C, 2.0)

# Plot the resulting PDF
using Plots
plot(x->pdf(C, x), 1, 5, yscale=:log10, ylim=(1e-3, 1))
```

- The result is a new numerically specified PDF on an extended grid, wrapped as a `NumericallyIntegrable` distribution.
- Both input PDFs must be normalized (integrate to 1) and have finite support.
- The convolution uses FFTW for efficiency and supports both vector and distribution types.

See the docstring for `fft_convolve` for more details and options.

## Implementation Details

The package uses numerical integration (via `QuadGK.jl`) to normalize the PDF and compute the CDF. For sampling, it uses a binned approximation of the CDF inversion method with a linear interpolation scheme.

For distributions with infinite support, a tangent transformation is used to map the infinite interval to (-1,1) before binning.

The `Interpolated` constructor leverages the `Interpolations.jl` package to create distributions from user-provided functions and grid points. It supports both `Constant()` and `Linear()` interpolation degrees, allowing for different levels of smoothness in the resulting distribution.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the general feature implementation strategy and checklist.
