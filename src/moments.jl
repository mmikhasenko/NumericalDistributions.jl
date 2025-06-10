"""
    numerical_moment(model::NumericallyIntegrable, n::Int, μ = mean(model))

Core implementation for numerical computation of distribution moments using integration.
Uses the package's integration facilities to directly compute moments from the PDF
rather than relying on analytical formulas, allowing moments computation for arbitrary
numerical distributions.
"""
function numerical_moment(model::NumericallyIntegrable, n::Int, μ = mean(model))
    return integral(x -> (x - μ)^n * pdf(model, x), model.support...)
end

"""
    mean(model::NumericallyIntegrable)

Numerical computation of distribution mean through direct integration.
Implements the Distributions.jl interface by computing the first raw moment
(using 0 as the reference point).
"""
mean(model::NumericallyIntegrable) = numerical_moment(model, 1, 0)

"""
    var(model::NumericallyIntegrable)

Numerical calculation of distribution variance using integration.
Implements the Distributions.jl interface by directly computing the second central moment
relative to the distribution's mean.
"""
var(model::NumericallyIntegrable) = numerical_moment(model, 2)

"""
    skewness(model::NumericallyIntegrable)

Implementation of skewness calculation for numerical distributions.
Computes the standardized third central moment by integrating over the PDF
and dividing by the cube of the standard deviation for proper normalization.
"""
function skewness(model::NumericallyIntegrable)
    σ = std(model)
    μ3 = numerical_moment(model, 3)
    return μ3 / σ^3
end

"""
    kurtosis(model::NumericallyIntegrable)

Numerical excess kurtosis implementation for arbitrary distributions.
Follows Distributions.jl convention by computing the normalized fourth central moment and subtracting 3 to make the normal distribution's kurtosis equal to zero.
"""
function kurtosis(model::NumericallyIntegrable)
    σ = std(model)
    μ4 = numerical_moment(model, 4)
    return μ4 / σ^4 - 3
end
