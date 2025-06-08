using NumericalDistributions
using NumericalDistributions.Interpolations

f(x) = (1 - x^2) * exp(-x / 2)
const n_sampling_bins = 1000
d1 = NumericallyIntegrable(f, (-0.5, 0.9); n_sampling_bins)
d2 = Interpolated(f, range(-0.5, 0.9, n_sampling_bins + 1); degree = Constant())
d3 = Interpolated(f, range(-0.5, 0.9, n_sampling_bins + 1); degree = Linear())

using BenchmarkTools

@btime rand($d1, 1000); # about 306.458 μs
@btime rand($d2, 1000); # about 2.029 s

# clear why slower
@btime rand($d3, 1000); # about 690.124 ms
