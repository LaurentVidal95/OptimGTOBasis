# Beware some imaginary parts are non negligeable.
function reference_eigenvectors(data::Dict{String, Any})
    ΨA = data["orba.re"] .+ im .* data["orba.im"]
    ΨB = data["orbb.re"] .+ im .* data["orbb.im"]
    ΨA, ΨB
end

function overlap(grid::QuadGrid{T1}, AOs_on_grid::Matrix{T2}) where {T1, T2<:Real}
    dot(grid, AOs_on_grid, AOs_on_grid)
end

######################### TODO: Code criterion j to optimize.

"""
For now just handles L2 projection
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
"""
function j_diatomic(A::Element{T1}, B::Element{T1},
                    RA::Vector{T2},  RB::Vector{T2},
                    Ψ_ref::Matrix{T3}, grid::QuadGrid{T2}) where {T1,T2 <: Real, T3}
    # Construct the AO_basis and eval on the integration grid
    AOs = vcat(construct_AOs(A; position=RA, grid.mmax, verbose=false),
               construct_AOs(B; position=RB, grid.mmax, verbose=false))
    𝐗 = eval_AOs(grid, AOs) ####################################### <- speedup needed
    
    # Compute the projection of ΨA on the AO basis
    n = length(𝐗)
    S = overlap(grid, 𝐗)
    # Switch from Dual type to Float64 if needed for inversion of S
    if eltype(𝐗) ≠ T2
        S = [x.value for x in overlap(grid, 𝐗)]
    end
    Γ = dot(grid, 𝐗, Ψ_ref)        
    C = inv(Symmetric(S))*Γ # projection coefficients

    sum(real([1 - 2*a'b + a'S*a for (a,b) in  zip(eachcol(C), eachcol(Γ))])) # sum all distances
end
