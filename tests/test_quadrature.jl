"""
Assert that the integration error over the given grid is lower
 that machine (double) precision
"""
function test_grid(grid::QuadGrid, ζ)
    g(X) = (π/ζ)^(-3/2)*exp(-ζ*norm(X)^2)
    ∫g = sum(grid.weights .* g.(grid.points))
    return abs(1 - ∫g) < 1e-15
end
