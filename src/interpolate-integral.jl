
function integral(
    itr::Interpolations.GriddedInterpolation{T,N,TC,Gridded{ST},K},
    a::Real,
    b::Real,
) where {T,N,TC,ST<:Interpolations.Constant,K}
    a > b && return -integral(itr, b, a)
    #
    grid = itr.knots[1]
    values = itr.coefs
    n = length(grid)

    # Handle out-of-bounds
    @boundscheck ((grid[1] <= a <= grid[end]) || Base.throw_boundserror(itr, a))
    @boundscheck ((grid[1] <= b <= grid[end]) || Base.throw_boundserror(itr, b))

    # Construct bin edges - these are not equally spaced
    edges = Vector{eltype(grid)}(undef, n + 1)
    edges[1] = grid[1]
    for i = 2:n
        edges[i] = 0.5 * (grid[i] + grid[i-1])
    end
    edges[n+1] = grid[n]

    idx_a = searchsortedlast(edges, a)
    idx_b = searchsortedlast(edges, b)
    # when the number is at the border, use the previous bin
    idx_b += (idx_b == n + 1) ? -1 : 0
    #
    s = zero(a)
    # integral of a is fully taken,
    # then subtract the part that should not be included
    # for the interval b, integrate from edge to b
    for i = idx_a:idx_b-1
        s += (edges[i+1] - edges[i]) * values[i]
    end
    s -= (a - edges[idx_a]) * values[idx_a]
    s += (b - edges[idx_b]) * values[idx_b]
    return s
end

function integral(
    itr::Interpolations.GriddedInterpolation{T,N,TC,Gridded{ST},K},
    a::Real,
    b::Real,
) where {T,N,TC,ST<:Interpolations.Linear,K}
    a > b && return -integral(itr, b, a)
    a == b && return zero(a)  # Early return for zero range integration

    grid = itr.knots[1]
    values = itr.coefs
    n = length(grid)

    # Handle out-of-bounds
    @boundscheck ((grid[1] <= a <= grid[end]) || Base.throw_boundserror(itr, a))
    @boundscheck ((grid[1] <= b <= grid[end]) || Base.throw_boundserror(itr, b))

    # Find indices of a and b in the grid
    idx_a = searchsortedlast(grid, a)
    idx_b = searchsortedlast(grid, b)
    idx_b += (idx_b == n) ? -1 : 0

    # Calculate result using full bins and adjusting for partial bins
    s = zero(promote_type(typeof(a), typeof(b)))
    # integral of a is fully taken,
    # then subtract the part that should not be included
    # for the interval b, integrate from edge to b
    for i = idx_a:idx_b-1
        s += in_bin_linear_integrate(grid[i], grid[i+1], values[i], values[i+1], grid[i+1])
    end
    s -= in_bin_linear_integrate(
        grid[idx_a],
        grid[idx_a+1],
        values[idx_a],
        values[idx_a+1],
        a,
    )
    s += in_bin_linear_integrate(
        grid[idx_b],
        grid[idx_b+1],
        values[idx_b],
        values[idx_b+1],
        b,
    )
    return s
end

function in_bin_linear_integrate(x1, x2, y1, y2, x)
    h = x2 - x1
    m = (y2 - y1) / h
    return (x - x1) * y1 + 0.5 * m * (x - x1)^2
end
