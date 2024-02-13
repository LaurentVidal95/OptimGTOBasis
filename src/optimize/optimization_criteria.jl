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
    norm_type             ::Symbol  # norm for the projection (L² or H¹)
    grids                 ::Vector{QuadGrid}
    gridtol               ::Float64
    interatomic_distances ::AbstractVector
end
function ProjectionCriterion(ref_data; gridtol=1e-7, norm_type=:L²)
    @assert norm_type ∈ (:L²,:H¹) # Only L² or H¹ norms
    ProjectionCriterion(ref_data.Ψs_ref, ref_data.TΨs_ref, norm_type,
                        ref_data.grids,  gridtol, (ref_data.Rhs .*2))
end

function Base.show(io::IO, crit::ProjectionCriterion)
    println(io, "Projection criterion for the $(crit.norm_type) norm")
end

function objective_function(criterion::ProjectionCriterion, A₀::Element,
                            B₀::Element, X::T...) where {T<:Real}
    Y = collect(X)
    if norm(Y) > 1e3
        foo = eltype(Y)==Float64 ? Y : map(x->x.value, Y)
        @info "High trial point norm: $(foo)"
    end
    J_Rhs = map(enumerate(criterion.interatomic_distances)) do (i, R)
        j_L2_diatomic(Y, A₀, B₀, R,
                      criterion.reference_functions[i], criterion.reference_kinetics[i],
                      criterion.grids[i]; criterion.norm_type)
    end
    sum(J_Rhs) / length(criterion.interatomic_distances)
end

"""
For now just handles L2 projection.
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
"""
function j_L2_diatomic(A::Element{T1}, B::Element{T1}, R::T2,
                       Ψ_ref::Matrix{T3}, TΨ_ref::Matrix{T3},
                       grid::QuadGrid{T2}; norm_type=:L²) where {T1,T2 <: Real, T3}

    # Construct the AO_basis and eval on the integration grid
    RA = [0.0, 0.0, -R/2];  RB = [0.0, 0.0, R/2]
    AOs = vcat(AO_basis(A; position=RA, grid.mmax, verbose=false),
               AO_basis(B; position=RB, grid.mmax, verbose=false))

    # Compute L² or H¹ projection on the AO basis
    C = eval_AOs(grid, AOs)
    M = norm_type==:L² ? dot(grid, C, C) : H¹_overlap(A, R)
    # Sanity check on the overlap
    if cond(M) > 1e5
        foo = cond(M)
        bar = typeof(foo)==Float64 ? foo : foo.value
        @warn "Overlap conditioning: $(bar)"
    end
    Mm12 = inv(sqrt(Symmetric(M)))
    C⁰ = C*Mm12 # AOs on the grid in orthonormal convention

    # Compute projection
    Π = dot(grid, C⁰, Ψ_ref)
    (norm_type==:H¹) && (Π .+= 2*dot(grid, C⁰, TΨ_ref))
    sum(1 .- Π'Π)
end
function j_L2_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                       R::T2, Ψ_ref::Matrix{T3}, TΨ_ref::Matrix{T3},
                       grid::QuadGrid{T2};
                       norm_type=:L²
                       ) where {T1,T2 <: Real, T3}
    nA = length(vec(A₀))
    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    (norm_type==:H¹) && (@assert A₀==B₀)

    XA, XB = (length(X)==nA) ? (X, X) : (X[1:nA], X[nA+1:end])
    A = Element(XA, A₀)
    B = Element(XB, B₀)
    # return j to minimize
    j_L2_diatomic(A, B, R, Ψ_ref, TΨ_ref, grid; norm_type)
end

struct EnergyCriterion{T<:Real} <: OptimizationCriterion
    reference_energies::Vector{T}
    interatomic_distances::Vector{T}
end
EnergyCriterion(ref_data; kwargs...) = EnergyCriterion(ref_data.Energies, (ref_data.Rhs .*2))

function objective_function(criterion::EnergyCriterion, A₀::Element, B₀::Element, X::T...) where {T<:Real}
    Y = collect(X)
    # Compute criterion
    J_Rhs = map(zip(criterion.reference_energies, criterion.interatomic_distances)) do (E, Rh)
        j_E_diatomic(Y, A₀, B₀, Rh/2, E)
    end
    sum(J_Rhs) / length(criterion.reference_energies)
end
function grad_objective_function!(criterion::EnergyCriterion, A₀::Element, B₀::Element, ∇J, X::T...) where {T<:Real}
    Y = collect(X)
    ∇Y = zero(Y)
    for (E, R) in zip(criterion.reference_energies, criterion.interatomic_distances)
            ∇Y .+= ∇j_E_diatomic(Y, A₀, B₀, R/2, E)
        end
    for i in 1:length(Y)
        ∇J[i] = ∇Y[i]
    end
    return ∇Y
end
function grad_objective_function(criterion::EnergyCriterion, A₀::Element, B₀::Element, X::T...) where {T<:Real}
    Y = collect(X)
    ∇Y = zero(Y)
    grad_objective_function!(criterion, A₀, B₀, X...)
    return ∇Y
end

function j_E_diatomic(A::Element{T1}, B::Element{T1}, Rh::T2,
    e_ref::T3) where {T1,T2,T3 <: Real}
    atom = "$(A.name) $(zero(Rh)) $(zero(Rh)) $(-Rh);
            $(B.name) $(zero(Rh)) $(zero(Rh)) $(Rh)"
    basis = A==B ? basis_string([A]) : basis_string([A,B])
    # Run pyscf
    mol = pyscf.M(;atom, basis,
                  symmetry=false,
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
function ∇j_E_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
                       Rh::T2, e_ref::T3) where {T1,T2,T3 <: Real}
    ∇j = []
    len_X = length(X)
    h = 1e-4
    for i in 1:len_X
        Xᵢph = X .+ [zeros(i-1)..., (h/2), zeros(len_X-i)...]
        Xᵢmh = X .- [zeros(i-1)..., (h/2), zeros(len_X-i)...]
        ∂i_j = (1/h)*(j_E_diatomic(Xᵢph, A₀, B₀, Rh, e_ref) -
                      j_E_diatomic(Xᵢmh, A₀, B₀, Rh, e_ref))
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
