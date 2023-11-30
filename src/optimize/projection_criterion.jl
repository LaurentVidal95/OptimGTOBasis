using GenericLinearAlgebra

"""
For now just handles L2 projection.
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
"""
function j_L2_diatomic(A::Element{T1}, B::Element{T1},
                       RA::Vector{T2},  RB::Vector{T2},
                       Ψ_ref::Matrix{T3}, grid::QuadGrid{T2}) where {T1,T2 <: Real, T3}
    # Construct the AO_basis and eval on the integration grid
    AOs = vcat(AO_basis(A; position=RA, grid.mmax, verbose=false),
               AO_basis(B; position=RB, grid.mmax, verbose=false))
    normalize_col(tab) = hcat(normalize.(eachcol(tab))...)
    # TODO: speedup
    𝐗 = normalize_col(eval_AOs(grid, AOs))

    # Compute the projection of ΨA on the AO basis
    S = dot(grid, 𝐗, 𝐗)
    # # Debug
    if eltype(𝐗) ≠ T2
    #     @show cond(S).value
    #     @show criterion.value
        S = [x.value for x in S]
    end
    Γ = dot(grid, 𝐗, Ψ_ref)
    C = inv(Symmetric(S))*Γ  # projection coefficients

    # Return sum of distances
    criterion = sum(norm(Ψ_ref_i - 𝐗*Ci)^2 for (Ψ_ref_i, Ci) in zip(eachcol(Ψ_ref), eachcol(C)))
    criterion # + 1e-5*cond(S)
end
function j_L2_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                       RA::Vector{T2},  RB::Vector{T2},
                       Ψ_ref::Matrix{T3}, grid::QuadGrid{T2},
                       ) where {T1,T2 <: Real, T3}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    XA, XB = (length(X)==nA) ? (X, X) : (X[1:nA], X[nA+1:end])
    A = Element(XA, A₀)
    B = Element(XB, B₀)
    # return j to minimize
    j_L2_diatomic(A, B, RA, RB, Ψ_ref, grid)
end

function j_E_diatomic(A::Element{T1}, B::Element{T1},
                      RA::Vector{T2},  RB::Vector{T2},
                      e_tot_ref::T3) where {T1,T2,T3 <: Real}
    atom = "$(A.name) $(RA[1]) $(RA[2]) $(RA[3]);
            $(B.name) $(RB[1]) $(RB[2]) $(RB[3])"
    basis = A==B ? basis_string([A]) : basis_string([A,B])
    # Run pyscf
    mol = pyscf.M(;atom, basis,
                  symmetry=false,
                  unit="bohr",
                  verbose=0,
                  )
    rhf_ground_state = mol.RHF().run()
    e_tot = rhf_ground_state.e_tot
    abs(e_tot - e_tot_ref)^2
end
function j_E_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                      RA::Vector{T2},  RB::Vector{T2},
                      e_tot_ref::T3) where {T1,T2,T3 <: Real}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    XA, XB = (length(X)==nA) ? (X, X) : (X[1:nA], X[nA+1:end])
    A = Element(XA, A₀)
    B = (length(X)==nA) ? A : Element(XB, B₀)

    # return j to minimize
    j_E_diatomic(A, B, RA, RB, e_tot_ref)
end
function ∇j_E_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                       RA::Vector{T2},  RB::Vector{T2},
                       e_tot_ref::T3) where {T1,T2,T3 <: Real}
    ∇j = []
    len_X = length(X)
    h = 1e-4
    for i in 1:len_X
        Xᵢph = X .+ [zeros(i-1)..., (h/2), zeros(len_X-i)...]
        Xᵢmh = X .- [zeros(i-1)..., (h/2), zeros(len_X-i)...]
        ∂i_j = (1/h)*(j_E_diatomic(Xᵢph, A₀, B₀, RA, RB, e_tot_ref) -
                      j_E_diatomic(Xᵢmh, A₀, B₀, RA, RB, e_tot_ref))
        push!(∇j, ∂i_j)
    end
    ∇j
end

# function legendre_coeffs(ζs, N_prim)
#     function f2opt(βs::Vector{T}) where {T<:Real}
#         ζs_approx = [_legendre_polynomial_to_coeff(i, N_prim, βs) for i in 1:N_prim]
#         norm(ζs .- ζs_approx)
#     end
#     optimize(f2opt, rand(length(ζs)), LBFGS())
# end

# function _legendre_polynomial_to_coeff(j::TI, N_prim::TI, βs::Vector{TR}) where {TI<:Int, TR<:Real}
#     exp(sum([β*Pl(int_to_reduced(j, 1, N_prim), iβ) for (iβ, β) in enumerate(βs)]))
# end

# function int_to_reduced(x, a, b)
#     @assert a ≤ x ≤ b "x ∉ [$a,$b]"
#     2*(x-a)/(b-a) - 1
# end
