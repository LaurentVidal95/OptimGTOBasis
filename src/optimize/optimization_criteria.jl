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

function orthogonal_projection(A::Element{T1}, B::Element{T1}, R::T2,
                               Ψ::Matrix{T3}, TΨ::Matrix{T3},
                               grid::QuadGrid{T2}; norm_type=:L²) where {T1,T2 <: Real, T3}
    C = eval_AOs(grid, A, B, R)
    M = overlap(grid, A, R; norm_type)
    Mm12 = inv(sqrt(Symmetric(M)))
    C⁰ = C*Mm12

    # Orthogonal projection for the given norm
    N_ao = size(C,2)
    Π = zeros(N_ao, N_ao)
    begin
        (norm_type==:L²) && (Π = dot(grid, C⁰, Ψ))
        (norm_type==:H¹) && (Π = dot(grid, C⁰, Ψ) + 2*dot(grid, C⁰, TΨ))
    end
    C⁰*Π
end

"""
Since ΨA and ΨB are equal for the current test casses I only put
Ψ_ref instead of ΨA, ΨB as argument.
T1 is the ForwardDiff compatible type, T2 is to be Float64 and T3 may be complex.
"""
function j_proj_diatomic(A::Element{T1}, B::Element{T1}, R::T2,
                         Ψ::Matrix{T3}, TΨ::Matrix{T3},
                         grid::QuadGrid{T2}; norm_type=:L²) where {T1,T2 <: Real, T3}
    # Compute L² or H¹ projection on the AO basis
    C = eval_AOs(grid, A, B, R)
    M = overlap(grid, A, R; norm_type)

    # Ψ_s_norm
    Ψ_norm = dot(grid, Ψ, Ψ)
    (norm_type==:H¹) && (Ψ_norm += 2*dot(grid, Ψ, TΨ))
    # Sanity check on the overlap
    if cond(M) > 1e5
        foo = eltype(M)==Float64 ? vec(A) : map(x->x.value, vec(A))
        bar = eltype(M)==Float64 ? M : map(x->x.value, M)
        @show foo
        @show bar
        @warn "Overlap conditioning: $(cond(foo))"
    end
    Mm12 = inv(sqrt(Symmetric(M)))
    C⁰ = C*Mm12

    # Compute projection
    Π = zero(Ψ_norm)
    begin
        (norm_type==:L²) && (Π = dot(grid, C⁰, Ψ))
        (norm_type==:H¹) && (Π = dot(grid, C⁰, Ψ) + 2*dot(grid, C⁰, TΨ))
    end
    sum(diag(Ψ_norm .- sum(Π'Π))) / length(Ψ)
end
function j_proj_diatomic(X::Vector{T1}, A₀::Element{T2}, B₀::Element{T2},
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
    j_proj_diatomic(A, B, R, Ψ_ref, TΨ_ref, grid; norm_type)
end

struct EnergyCriterion{T<:Real} <: OptimizationCriterion
    reference_energies::Vector{T}
    interatomic_distances::Vector{T}
end
EnergyCriterion(ref_data; kwargs...) =
    EnergyCriterion(ref_data.Energies, (ref_data.Rhs .*2))

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


function objective_function(criterion::ProjectionCriterion, A₀::Element,
                            B₀::Element, X::T...) where {T<:Real}
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
                      criterion.norm_type)
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
