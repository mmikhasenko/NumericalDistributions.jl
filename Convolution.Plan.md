# Convolution.Plan.md

## Status Summary (June 2025)
- **Core FFT-based convolution feature is implemented, tested, and integrated.**
- **All major steps are complete.**
- Optional/documentation/edge-case improvements remain.

---

## Numerical Convolution via FFTW: Implementation Plan

### 1. Overview
Implement a function to numerically convolve two functions defined on the same grid (with equal step size) using FFTW. The result should be a new function on an extended grid, preserving the step size.

### 2. Incorporate Prototype Logic (with Interpolated)
- [x] Use the working algorithm from `convolution_prototype.jl` as the reference implementation.
- [x] Instead of the prototype's `OnGrid`, use the package's `Interpolated` constructor for tabulated/interpolated PDFs.
- [x] The prototype provides:
  - `pdf_convolve(f::AbstractVector, g::AbstractVector; Δ, t0_f, t0_g, pow2)`
  - `pdf_convolve(f::OnGrid, g::OnGrid; pow2)`
- [x] In the package, adapt these to:
  - `convolve_pdf(pdf1::AbstractVector, pdf2::AbstractVector; Δ, t0_1, t0_2, pow2)`
  - `convolve_pdf(pdf1::Interpolated, pdf2::Interpolated; pow2, degree)`
- [x] Ensure the result is an `Interpolated` object, using the same grid logic and interpolation degree as the inputs.
- [x] Add Julia docstrings and type annotations.
- [x] The docstring must clearly state that both `pdf1` and `pdf2` must be normalized probability density functions.

### 3. Implementation Steps

#### 3.1. Add FFTW as a dependency [COMPLETED]
- [x] Update `Project.toml` to include `FFTW` in `[deps]` and `[compat]`.

#### 3.2. Implement convolution function [COMPLETED]
- [x] Create a new file `convolution.jl` in `src/`.
- [x] Implement the convolution functions based on the prototype, but using `Interpolated` instead of `OnGrid`:
  - `convolve_pdf(pdf1::AbstractVector, pdf2::AbstractVector; Δ, t0_1, t0_2, pow2)`
  - `convolve_pdf(pdf1::Interpolated, pdf2::Interpolated; pow2, degree)`
- [x] Ensure the result is an `Interpolated` object for further use.
- [x] The docstring must warn that both arguments must be normalized PDFs.

#### 3.3. Expose convolution in the main module [COMPLETED]
- [x] Add `include("convolution.jl")` and export the convolution function in `NumericalDistributions.jl`.

#### 3.4. Write tests [COMPLETED]
- [x] Add tests in `test/runtests.jl`:
  - [x] Test convolution of simple functions (e.g., box, Gaussian) using both vector and `Interpolated` APIs.
  - [x] Compare with analytical results where possible.
  - [x] Test edge cases (different lengths, zeros, etc.).
  - [x] Adjust tolerances and fix grid alignment as needed.

### 4. Documentation [COMPLETED]
- [x] Add docstrings to all new functions.
- [x] Briefly document the feature in the README.

### 5. (Optional) API Extensions [FUTURE]
- [ ] Support for convolution of distributions or interpolated objects with different grid types.
- [ ] Support for non-uniform grids (future work).

---

**Next Steps:**
- [x] Add a usage example to the README.
- [x] Add more edge-case or mixed-type tests (including with Distributions.jl types).
- [ ] (Optional) Further generalize API for future needs.
- [x] Add error/warning for infinite support and large grid sizes in convolution.
