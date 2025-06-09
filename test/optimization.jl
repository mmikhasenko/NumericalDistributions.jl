using NumericalDistributions
using NumericalDistributions.Interpolations

f(x) = (1 - x^2) * exp(-x / 2)
const n_sampling_bins = 1000
d1 = NumericallyIntegrable(f, (-0.5, 0.9); n_sampling_bins)
d2 = interpolated(f, range(-0.5, 0.9, n_sampling_bins + 1); degree = Constant())
d3 = interpolated(f, range(-0.5, 0.9, n_sampling_bins + 1); degree = Linear())

using BenchmarkTools

@btime rand($d1, 1000); # 34.542 μs
@btime rand($d2, 1000); # 28.500 μs

# clear why slower
@btime rand($d3, 1000); # 754.542 μs
