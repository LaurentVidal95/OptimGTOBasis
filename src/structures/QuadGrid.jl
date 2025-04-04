import LinearAlgebra.dot

"""
Numerical integration grid
"""
struct QuadGrid{T<:Real}
    points  ::AbstractArray{Vector{T}} # vector of cartesian coordinates in 𝐑³
    weights ::Vector{T}                # corresponding weights
    # maximum quantum number m that can be integrated
    # Integrating spherical harmonics with |m| > mmax
    # might result in uncontrolled errors.
    mmax    ::Int64
end

"""
Construct the integration grid from a data dictionary.
"""
function QuadGrid(data::Dict{String, Any})
    # Construct grid
    Rh = data["Rh"]
    mmax = Int64(data["mmax"])
    function polsph_to_cart(μ, cos_ν, φ)
        common = Rh*sinh(μ)*√(1-cos_ν^2)
        [common*cos(φ), common*sin(φ), Rh*cosh(μ)*cos_ν]
    end
    weights = vec(data["dV"])
    points = [polsph_to_cart(μ, cos_ν, φ) for (μ, cos_ν, φ) in
              zip(vec(data["mu"]), vec(data["cth"]), vec(data["phi"]))]
    QuadGrid(points, weights, mmax)
end

"""
Dot product on the numerical quadrature grid defined as:

    ⟨ψ1|ψ2⟩ = ∑_k ω_k * ψ1_k * ψ2_k

where [ψ1]_k and [ψ2]_k are the vectors of two functions ψ1 and ψ2
evaluated on the quadrature points.
"""
function dot(grid::QuadGrid{T1}, Ψ1, Ψ2) where {T1<:Real}
    (grid.weights .* Ψ1)'Ψ2
end

"""
Assert that the integration error over the given grid is lower
than machine (double) precision
"""
function test_spread(grid::QuadGrid{T}, ζ, Rh::T; gridtol) where {T<:Real}
    g(X) = (π/ζ)^(-3/2)*exp(-ζ*norm(X .- [0,0,Rh])^2)
    ∫g = sum(grid.weights .* g.(grid.points))
    return abs(1 - ∫g) < gridtol
end

"""
Computes within 5 units the maximum ζ that can be integrated on the integration grid
with maximum error gridtol.
"""
function compute_spread_lim(grid::QuadGrid{T}, Rh::T; gridtol) where {T<:Real}
    spread_lim = 20
    while test_spread(grid, spread_lim, Rh; gridtol)
        spread_lim += 5
    end
    spread_lim
end
