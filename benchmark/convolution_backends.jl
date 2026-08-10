#!/usr/bin/env julia
# Run from package root: julia --project=. benchmark/convolution_backends.jl

using NumericalDistributions
using Distributions
using FFTW
using BenchmarkTools
using Statistics
using Printf
using ReverseDiff

function make_problem(n1::Int, n2::Int)
    Δ = 0.01
    x1 = range(-2, 4; length = n1)
    x2 = range(-3, 3; length = n2)
    y1 = pdf.(Uniform(-0.5, 3.5), x1)
    y2 = pdf.(Normal(0, 0.3), x2)
    return y1, y2, Δ, first(x1), first(x2)
end

function max_rel_err(a, b)
    denom = max(maximum(abs.(b)), eps())
    return maximum(abs.(a .- b)) / denom
end

function accuracy_table()
    println("\n## Accuracy (max relative error vs direct convolution)\n")
    println("| n1 | n2 | :fftw vs direct | :generic vs direct | :fftw vs :generic |")
    println("|----|----|-----------------|--------------------|-------------------|")
    for (n1, n2) in ((50, 50), (200, 200), (500, 500), (1000, 1000))
        y1, y2, Δ, t0_1, t0_2 = make_problem(n1, n2)
        kw = (; Δ, t0_1, t0_2)
        _, h_dir = NumericalDistributions._direct_convolve(y1, y2; kw...)
        _, h_fftw = NumericalDistributions._fft_convolve(y1, y2; kw...)
        _, h_gen = NumericalDistributions._generic_fft_convolve(y1, y2; kw...)
        e_fd = max_rel_err(h_fftw, h_dir)
        e_gd = max_rel_err(h_gen, h_dir)
        e_fg = max_rel_err(h_fftw, h_gen)
        @printf("| %d | %d | %.3e | %.3e | %.3e |\n", n1, n2, e_fd, e_gd, e_fg)
    end
end

function timing_table()
    println("\n## Timing (median of BenchmarkTools, 1 sample + 10 evals)\n")
    println("| n1 | n2 | :fftw (μs) | :generic (μs) | :direct (μs) | generic/fftw | direct/fftw |")
    println("|----|----|------------|---------------|--------------|--------------|-------------|")
    for (n1, n2) in ((100, 100), (300, 300), (500, 500), (1000, 1000))
        y1, y2, Δ, t0_1, t0_2 = make_problem(n1, n2)
        kw = (; Δ, t0_1, t0_2)
        t_fftw = @belapsed NumericalDistributions._fft_convolve($y1, $y2; $kw...)
        t_gen = @belapsed NumericalDistributions._generic_fft_convolve($y1, $y2; $kw...)
        t_dir = @belapsed NumericalDistributions._direct_convolve($y1, $y2; $kw...)
        @printf(
            "| %d | %d | %.1f | %.1f | %.1f | %.1f× | %.1f× |\n",
            n1,
            n2,
            t_fftw * 1e6,
            t_gen * 1e6,
            t_dir * 1e6,
            t_gen / t_fftw,
            t_dir / t_fftw,
        )
    end
end

function tracked_timing()
    println("\n## ReverseDiff tracked vector (n1=n2=201, :auto → :direct)\n")
    y1, y2, Δ, t0_1, t0_2 = make_problem(201, 201)
    y2t = ReverseDiff.track.(y2)
    kw = (; Δ, t0_1, t0_2)
    t_auto = @belapsed NumericalDistributions._convolve_vectors($y1, $y2t; $kw...)
    t_dir = @belapsed NumericalDistributions._direct_convolve($y1, $y2t; $kw...)
    @printf(":auto:   %.1f μs\n", t_auto * 1e6)
    @printf(":direct: %.1f μs\n", t_dir * 1e6)
end

println("# Convolution backend benchmark")
println("Julia ", VERSION)
println("Host: ", Sys.MACHINE)
accuracy_table()
timing_table()
tracked_timing()
