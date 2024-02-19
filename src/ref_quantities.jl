function orthogonal_projection(A::Element, B::Element, R::T1,
                               Ψ::Matrix{T2}, TΨ::Matrix{T2},
                               grid::QuadGrid{T1}; norm_type=:L²) where {T1 <: Real, T2}
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

function quadrupole_moment(grid::QuadGrid, A::Element, B::Element,
                           Ψs::Matrix{T}, R::T;
                           verbose=true) where {T<:Real}
    # DEBUG: only for closed shell... Otherwise you have to track the alpha and beta part..
    (verbose) && (@warn "For closed shell only!")
    ρ = 2*sum(map(x->x .^2, eachcol(Ψs))) # density associated to Ψ

    quadmoment = zeros(3,3)
    # nuclear contribution for the zz quadripolar moment
    quadmoment[end] = (A.charge + B.charge)*(R^2)/4

    for i in 1:3
        for j in 1:3
            δᵢⱼ = iszero(i-j) ? one(i) : zero(i)
            quadmoment_n_terms = map(enumerate(grid.points)) do (ir, r)
                (1/2)*grid.weights[ir]*( 3*r[i]*r[j] - δᵢⱼ*norm(r)^2 )*ρ[ir]
            end
            quadmoment[i, j] -= sum(quadmoment_n_terms)
        end
    end
    quadmoment
end

function dipole_moment()
    nothing #TODO
end

function polarisability(mol::PyObject; ε=1e-4, kwargs...)
    # Extract data from the mol object
    A = mol.atom_symbol(0) # element 1
    B = mol.atom_symbol(1) # element 2
    R = norm(mol.atom_coord(1) - mol.atom_coord(0))

    # Add running Helfem and computing finite difference of the energy with constant
    # electronic field in the z direction.
    α = zero(Int64)

    # TODO

    α
end
