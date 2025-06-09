function numerical_moment(model::NumericallyIntegrable, n::Int, μ = mean(model))
    return integral(x -> (x - μ)^n * pdf(model, x), model.support...)
end

mean(model::NumericallyIntegrable) = numerical_moment(model, 1, 0)

var(model::NumericallyIntegrable) = numerical_moment(model, 2)

function skewness(model::NumericallyIntegrable)
    σ = std(model)
    μ3 = numerical_moment(model, 3)
    return μ3 / σ^3
end

function kurtosis(model::NumericallyIntegrable)
    σ = std(model)
    μ4 = numerical_moment(model, 4)
    return μ4 / σ^4 - 3
end
