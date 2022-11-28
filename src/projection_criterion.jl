"""
For now just handles L2 projection
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
"""
function j_L2_diatomic(A::Element{T1}, B::Element{T1},
                       RA::Vector{T2},  RB::Vector{T2},
                       Ψ_ref::Matrix{T3}, grid::QuadGrid{T2}) where {T1,T2 <: Real, T3}
    # Construct the AO_basis and eval on the integration grid
    AOs = vcat(construct_AOs(A; position=RA, grid.mmax, verbose=false),
               construct_AOs(B; position=RB, grid.mmax, verbose=false))
    𝐗 = eval_AOs(grid, AOs) ####################################### <- speedup needed

    # Compute the projection of ΨA on the AO basis
    n = length(𝐗)
    S = dot(grid, 𝐗, 𝐗)

    # Switch from Dual type to Float64 if needed for inversion of S
    if eltype(𝐗) ≠ T2
        S = [x.value for x in S]
    end
    Γ = dot(grid, 𝐗, Ψ_ref)
    C = inv(Symmetric(S))*Γ # projection coefficients

    sum(real([1 - 2*a'b + a'S*a for (a,b) in  zip(eachcol(C), eachcol(Γ))])) # sum all distances
end
function j_L2_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                       RA::Vector{T2},  RB::Vector{T2},
                       Ψ_ref::Matrix{T3}, grid::QuadGrid{T2},
                       ) where {T1,T2 <: Real, T3}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    XA, XB = X[1:nA], X[nA+1:end]
    (eltype(X) ≠ T2) && (@show [x.value for x in X])
    A = Element(XA, A₀)
    B = Element(XB, B₀)
    # return j to minimize
    j_L2_diatomic(A, B, RA, RB, Ψ_ref, grid)
end
