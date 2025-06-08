### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ c073bfc3-e473-4c58-a7cf-15777c90d15c
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    # 
    using NumericalDistributions
    # using NumericalDistributions.Distributions
    using Plots
end

# ╔═╡ 7e78061f-e2b2-41be-bf60-a7b3acf1675d
theme(:boxed)

# ╔═╡ aa4d233a-9b00-495d-aa62-94b35d7385e6
begin
    # Define two custom PDFs on different supports
    x1 = -2:0.01:2
    x2 = -1:0.01:8
end

# ╔═╡ ee2d7f4d-ba1b-4bd7-aa56-5db35e8fcd40
f1(x) = exp(-x^4)                    # Even, super-Gaussian

# ╔═╡ cc2158ab-1481-4f20-a0dd-f656960aa4c9
f2(x) = 1 / ((x - 3)^2 + 1)          # Cauchy-like, centered at 3

# ╔═╡ 46220ab4-4afc-48f9-a7b5-1da01d2e2623
md"""
Create interpolated distributions
"""

# ╔═╡ f0a61d86-bf7b-4c89-9aa7-40a93828ffc3
A = Interpolated(f1, x1);

# ╔═╡ df16ccf0-cf23-4c45-bd39-7318e17bd483
B = Interpolated(f2, x2);

# ╔═╡ 87ba4a95-c6f2-4fd3-b34f-d134a996d687
C = fft_convolve(A, B);

# ╔═╡ ad903773-d706-46e4-9b9d-50d87d77161f
# Evaluate the resulting PDF at a point
pdf(C, 2.0)

# ╔═╡ ae035511-1a07-4cd0-8327-87dc46466917
plot(x -> pdf(C, x), -3, 11, yscale = :log10, ylim = (1e-3, 1))

# ╔═╡ Cell order:
# ╠═c073bfc3-e473-4c58-a7cf-15777c90d15c
# ╠═7e78061f-e2b2-41be-bf60-a7b3acf1675d
# ╠═aa4d233a-9b00-495d-aa62-94b35d7385e6
# ╠═ee2d7f4d-ba1b-4bd7-aa56-5db35e8fcd40
# ╠═cc2158ab-1481-4f20-a0dd-f656960aa4c9
# ╟─46220ab4-4afc-48f9-a7b5-1da01d2e2623
# ╠═f0a61d86-bf7b-4c89-9aa7-40a93828ffc3
# ╠═df16ccf0-cf23-4c45-bd39-7318e17bd483
# ╠═87ba4a95-c6f2-4fd3-b34f-d134a996d687
# ╠═ad903773-d706-46e4-9b9d-50d87d77161f
# ╠═ae035511-1a07-4cd0-8327-87dc46466917
