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
"""
function j_diatomic(A::Element{T1}, B::Element{T1},
                    RA::Vector{T2},  RB::Vector{T2},
                    Ψ_ref::Matrix{T3}, grid::QuadGrid{T2}) where {T1,T2 <: Real, T3}
    # Construct the AO_basis and eval on the integration grid
    AOs = vcat(construct_AOs(A; position=RA, grid.mmax, verbose=false),
               construct_AOs(B; position=RB, grid.mmax, verbose=false))
    𝐗 = eval_AOs(grid, AOs) ####################################### <- speedup needed
    
    # Orthonormalize
    S = [x.value for x in overlap(grid, 𝐗)]
    𝐗_ortho = 𝐗*inv(sqrt(Symmetric(S)))

    # Compute projection on the AO basis
    j_out = zero(T1) # not sure about type here for ForwardDiff
    for X in eachcol(𝐗) # run over all AOs
        j_out += sum(abs2, dot(grid, X, Ψ_ref)) # project A eigenvectors
    end
    j_out  # - 1e-4*log(cond(S)) # add this constraint on conditioning ?
end
