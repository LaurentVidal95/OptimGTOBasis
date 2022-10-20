import LinearAlgebra.dot

struct QuadGrid{T<:Real}
    points  ::AbstractArray{T}
    weigths ::Vector{T}
end

"""
Adapter les arguments en fonction des besoins
Il peut être mieux de prendre directement un tableau en entré plutôt que f
"""
function ∫(grid::QuadGrid{T}, f::F) where {T<:Real, F}
    sum(grid.weigths .* f.(grid.points))
end

"""
Dot product on the numerical quadrature grid.

    ⟨ψ1|ψ2⟩ = ∑_k ω_k * ψ1_k * ψ2_k

ψ1 and ψ2 are typically two functions evaluated on the quadrature points
and vectorized.
"""
function dot(grid::QuadGrid{T}, Ψ1::Vector{T}, Ψ2::Vector{T}) where {T<:Real}
    @assert( length(grid.weights) == length(Ψ1) == length(Ψ2) )
    sum(grid.weigths .* ψ1 .* Ψ2)
end
