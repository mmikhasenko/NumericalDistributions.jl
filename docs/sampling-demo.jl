### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 6ae57d44-447b-11f0-374a-0f43aadcc968
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    # 
    using NumericalDistributions
    using Plots
end

# ╔═╡ 3b8f966b-4acc-4064-bed4-3b90d7281366
md"""
# Numerical sampling

This notebook demonstrates how different algorighms for sampling work.
When interpolated linearly, or interpolated with a const.
"""

# ╔═╡ f7b5cae2-abea-47eb-b4c7-270800a2705a
theme(:boxed)

# ╔═╡ e5958f19-040b-4503-a8fb-2565b27d575e
f(x) = (1 - x^2) * exp(-(x-0.1)^2 / 0.2)

# ╔═╡ 5ef71ef7-91e4-4292-bfb0-53780fd9f37a
x_range = (-0.5, 0.8)

# ╔═╡ 9dec7c61-df4e-4043-a1d4-9181b64e9174
const n_sampling_bins = 3

# ╔═╡ fba9fde0-307e-43ce-8133-e101c71439e6
md"""
## Current
Current interpolation of the sampling algorithm uses right bin center
"""

# ╔═╡ 20af1c3f-4cb6-42b7-a1d7-e72df0e5b980
d = NumericallyIntegrable(f, x_range; n_sampling_bins);

# ╔═╡ 9c49d97d-fd50-4eb1-b12b-b9467a4bf860
let
    data = rand(d, 1_000_000)
    bins = range(x_range..., 100)
    scale = length(data) * (bins[2] - bins[1])
	# 
    plot(x -> pdf(d, x) * scale, x_range..., lab="true")
    stephist!(data; bins, lab="sampled")
end

# ╔═╡ 5a60e4f6-fd8b-4d91-bdc2-de88b905cd39
md"""
## Interpolated Constant
uses half bin for the first and last bin, the others with a bin center
"""

# ╔═╡ 55d99ee7-aceb-4139-b7e2-9907c421a8f5
id = Interpolated(
    f,
    range(x_range..., n_sampling_bins + 1);
    degree = NumericalDistributions.Interpolations.Constant(),
);

# ╔═╡ 77346e05-eee8-477a-9e7e-30e596acf71a
let
    data = rand(id, 1_000_000)
    bins = range(x_range..., 100)
    scale = length(data) * (bins[2] - bins[1])
	# 
    plot(x -> pdf(d, x) * scale, x_range..., lab="true")
	plot!(x -> pdf(id, x) * scale, x_range..., lab="constant-interp.")
	stephist!(data; bins, lab="sampled")
end

# ╔═╡ 0aef928f-90b5-4eb6-847f-891f7d900acd
md"""
## Inv cdf
look very reasonable
"""

# ╔═╡ 6772ec63-6370-4815-a79e-c914c6f39058
il = Interpolated(
    f,
    range(x_range..., n_sampling_bins + 1);
    degree = NumericalDistributions.Interpolations.Linear(),
);

# ╔═╡ 6aa1af8d-3436-4209-b0b1-ed02141919fc
begin
	plot(aspect_ratio=1)
	plot!(x->cdf(id, x), range(x_range..., 100))
	plot!(x->cdf(il, x), range(x_range..., 100))
end

# ╔═╡ 2a935941-9ed8-469e-9edf-bec9506e9a9b
begin
	plot(aspect_ratio=1)
	plot!(x->invcdf(id, x), range(0, 1, 100))
	plot!(x->invcdf(il, x), range(0, 1, 100))
end

# ╔═╡ 011cea46-070c-417e-b664-a335adf4f1ff


# ╔═╡ Cell order:
# ╟─3b8f966b-4acc-4064-bed4-3b90d7281366
# ╠═6ae57d44-447b-11f0-374a-0f43aadcc968
# ╠═f7b5cae2-abea-47eb-b4c7-270800a2705a
# ╠═e5958f19-040b-4503-a8fb-2565b27d575e
# ╠═5ef71ef7-91e4-4292-bfb0-53780fd9f37a
# ╠═9dec7c61-df4e-4043-a1d4-9181b64e9174
# ╟─fba9fde0-307e-43ce-8133-e101c71439e6
# ╠═20af1c3f-4cb6-42b7-a1d7-e72df0e5b980
# ╠═9c49d97d-fd50-4eb1-b12b-b9467a4bf860
# ╟─5a60e4f6-fd8b-4d91-bdc2-de88b905cd39
# ╠═77346e05-eee8-477a-9e7e-30e596acf71a
# ╠═55d99ee7-aceb-4139-b7e2-9907c421a8f5
# ╟─0aef928f-90b5-4eb6-847f-891f7d900acd
# ╠═6aa1af8d-3436-4209-b0b1-ed02141919fc
# ╠═2a935941-9ed8-469e-9edf-bec9506e9a9b
# ╠═6772ec63-6370-4815-a79e-c914c6f39058
# ╠═011cea46-070c-417e-b664-a335adf4f1ff
