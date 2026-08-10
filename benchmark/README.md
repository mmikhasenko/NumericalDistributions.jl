# Benchmarks

Convolution backend comparison (accuracy and timing):

```bash
julia --project=. -e '
  using Pkg
  Pkg.add(["BenchmarkTools", "ReverseDiff"])
  include("benchmark/convolution_backends.jl")
'
```

These packages are only needed to run the script locally; they are not package dependencies.
