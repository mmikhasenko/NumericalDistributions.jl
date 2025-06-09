using NumericalDistributions
using Interpolations
using Distributions
using Statistics
using Random
using QuadGK
using Test
using ForwardDiff

include("test-basic.jl")
include("test-sampling.jl")
include("test-moments.jl")

include("test-interpolated.jl")
include("test-quantile.jl")

include("test-convolution.jl")
