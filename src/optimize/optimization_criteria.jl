using GenericLinearAlgebra

import Base.show

abstract type OptimizationCriterion end

"""
Optimize the distance between approach density and a reference
density ρ_ref on a discretization grid.
"""
struct ProjectionCriterion <: OptimizationCriterion
    reference_functions   ::Vector{Matrix}
    reference_kinetics    ::Vector{Matrix}
    norm_type             ::Symbol # norm for the projection (L² or H¹)
    grids                 ::Vector{QuadGrid}
    gridtol               ::Float64
    interatomic_distances ::AbstractVector
end
function ProjectionCriterion(ref_data; gridtol=1e-7, norm_type=:L²)
    @assert norm_type ∈ (:L²,:H¹) # Only L² or H¹ norms
    ProjectionCriterion(ref_data.reference_MOs, ref_data.reference_∇MOs, norm_type,
                        ref_data.grids, gridtol, ref_data.interatomic_distances)
end

function (crit::ProjectionCriterion)(basis_A::BasisSet, basis_B::BasisSet, i::Int)
    A = basis_A.element; B = basis_B.element;
    X = vec(A); Y = vec(B)
    Z = X==Y ? X : vcat(X,Y)

    j_proj_diatomic(Z, A, B, crit.interatomic_distances[i],
                    crit.reference_functions[i], crit.reference_kinetics[i],
                    crit.grids[i];
                    crit.norm_type)
end

function Base.show(io::IO, crit::ProjectionCriterion)
    println(io, "Projection criterion for the $(crit.norm_type) norm")
end

"""
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
    """
function j_proj_diatomic(A::Element{T1}, B::Element{T1}, R::T2,
                         Ψ::Matrix{T3}, TΨ::Matrix{T3},
                         grid::QuadGrid{T2}; norm_type=:L²,
                         shift=1.) where {T1,T2 <: Real, T3}
    # Compute L² or H¹ projection on the AO basis
    C = eval_AOs(grid, A, B, R)
    SA = overlap(grid, A, B, R; norm_type)
    N_ao = size(C, 2)
    N_mo = size(Ψ, 2)

    # Sanity check on the overlap
    if cond(SA) > 1e5
        foo = eltype(SA)==Float64 ? vec(A) : map(x->x.value, vec(A))
        bar = eltype(SA)==Float64 ? SA : map(x->x.value, SA)
        @warn "Bad overlap conditioning ($(cond(foo))) for interatomic distance $R a.u."
    end
    SA_inv = inv(Symmetric(SA))
    # Compute projection
    
    Π = zeros(eltype(C), N_ao, N_mo)
    begin
        (norm_type==:L²) && (Π = dot(grid, C, shift .* Ψ))
        (norm_type==:H¹) && (Π = dot(grid, C, shift .* Ψ) -
                             dot(grid, eval_AOs(grid, A, B, R; deriv=2), Ψ))
        # (norm_type==:H¹) && (Π = dot(grid, C, shift .* Ψ) + 2*dot(grid, C, TΨ)) # OLD (keep for debug)
    end
    -2*tr(Π'*SA_inv*Π)
end
function j_proj_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                         R::T2, Ψ_ref::Matrix{T3}, TΨ_ref::Matrix{T3},
                         grid::QuadGrid{T2};
                         norm_type=:L²,
                         shift=1
                         ) where {T1,T2 <: Real, T3}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    (norm_type==:H¹) && (@assert A₀==B₀)

    XA, XB = (length(X)==nA) ? (X, X) : (X[1:nA], X[nA+1:end])
    A = Element(XA, A₀)
    B = Element(XB, B₀)
    # return j to minimize
    j_proj_diatomic(A, B, R, Ψ_ref, TΨ_ref, grid; norm_type, shift)
end

struct EnergyCriterion{T<:Real} <: OptimizationCriterion
    reference_energies::Vector{T}
    interatomic_distances::Vector{T}
end
EnergyCriterion(ref_data; kwargs...) =
    EnergyCriterion(ref_data.energies, ref_data.interatomic_distances)

function (crit::EnergyCriterion)(basis_A::BasisSet, basis_B::BasisSet, i::Int)
    A = basis_A.element; B = basis_B.element;
    X = vec(A); Y = vec(B)
    Z = X==Y ? X : vcat(X,Y)
    
    j_E_diatomic(Z, A, B, crit.interatomic_distances[i]/2, crit.reference_energies[i])
end

function j_E_diatomic(A::Element{T1}, B::Element{T1}, Rh::T2,
    e_ref::T3) where {T1,T2,T3 <: Real}
    atom = "$(A.name) $(zero(Rh)) $(zero(Rh)) $(-Rh);
            $(B.name) $(zero(Rh)) $(zero(Rh)) $(Rh)"
    basis = A==B ? basis_string([A]) : basis_string([A,B])
    # Run pyscf
    mol = pyscf.M(;atom, basis,
                  unit="bohr",
                  verbose=0,
                  )
    rhf_ground_state = mol.RHF().run()
    abs(rhf_ground_state.e_tot - e_ref)^2
end
function j_E_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                      Rh::T2, e_ref::T3) where {T1,T2,T3 <: Real}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    XA, XB = (length(X)==nA) ? (X, X) : (X[1:nA], X[nA+1:end])
    A = Element(XA, A₀)
    B = (length(X)==nA) ? A : Element(XB, B₀)

    # return j to minimize
    j_E_diatomic(A, B, Rh, e_ref)
end
# function ∇j_E_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
#                        Rh::T2, e_ref::T3) where {T1,T2,T3 <: Real}
#     ∇j = []
#     len_X = length(X)
#     h = 1e-4
#     for i in 1:len_X
#         Xᵢph = X .+ [zeros(i-1)..., (h/2), zeros(len_X-i)...]
#         Xᵢmh = X .- [zeros(i-1)..., (h/2), zeros(len_X-i)...]
#         ∂i_j = (1/h)*(j_E_diatomic(Xᵢph, A₀, B₀, Rh, e_ref) -
#                       j_E_diatomic(Xᵢmh, A₀, B₀, Rh, e_ref))
#         push!(∇j, ∂i_j)
#     end
#     ∇j
# end


function objective_function(criterion::ProjectionCriterion, A₀::Element,
                            B₀::Element, X::T...;
                            kwargs...) where {T<:Real}
    Y = collect(X)

    # log optimization to ensure that the exponents are positive
    n_exps = len_exps(A₀)
    Y[1:n_exps] .= exp.(Y[1:n_exps])

    if norm(Y) > 1e3
        foo = eltype(Y)==Float64 ? Y : map(x->x.value, Y)
        @info "High trial point norm: $(foo)"
    end
    J_Rhs = map(enumerate(criterion.interatomic_distances)) do (i, R)
        j_proj_diatomic(Y, A₀, B₀, R,
                        criterion.reference_functions[i], criterion.reference_kinetics[i],
                        criterion.grids[i];
                        criterion.norm_type,
                        kwargs...)
    end
    sum(J_Rhs) / length(criterion.interatomic_distances)
end

function objective_function(criterion::EnergyCriterion, A₀::Element, B₀::Element,
                            X::T...) where {T<:Real}
    Y = collect(X)
    n_exps = len_exps(A₀)
    Y[1:n_exps] .= exp.(Y[1:n_exps])

    # Compute criterion
    J_Rhs = map(zip(criterion.reference_energies, criterion.interatomic_distances)) do (E, R)
        j_E_diatomic(Y, A₀, B₀, R/2, E)
    end
    sum(J_Rhs) / length(criterion.reference_energies)
end
