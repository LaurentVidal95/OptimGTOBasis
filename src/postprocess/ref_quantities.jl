function ref_properties(ref_basis_name, refdir)
    # First compute the reference HELFEM quantities
    data_helfem = read_helfem_data(ref_basis_name, refdir)
    # Extract system's name from data
    el1 = data_helfem.elements[1].name
    el2 = data_helfem.elements[2].name
    system = el1==el2 ? el1*"₂" : el1*el2

    ref_data = (; energy=data_helfem.energies,
                dipole=data_helfem.dipoles,
                quadrupole=data_helfem.quadrupoles,
                forces=.- data_helfem.forces,
                data_helfem.interatomic_distances,
                system)
    ref_data
end

function orthogonal_projection(A::Element, B::Element, R::T1,
                               Ψ::Matrix{T2}, TΨ::Matrix{T2},
                               grid::QuadGrid{T1}; norm_type=:L²) where {T1 <: Real, T2}
    C = eval_AOs(grid, A, B, R)
    SA = overlap(grid, A, B, R; norm_type)
    SA_inv = inv(Symmetric(SA))

    # Orthogonal projection for the given norm
    N_ao = size(C,2)
    N_mo = size(Ψ,2)
    Π = zeros(eltype(C), N_ao, N_mo)
    begin
        (norm_type==:L²) && (Π = dot(grid, C⁰, Ψ))
        (norm_type==:H¹) && (Π = dot(grid, C, shift .* Ψ) -
                             dot(grid, eval_AOs(grid, A, B, R; deriv=2), Ψ))
    end
    C*SA_inv*Π
end

function dipole_moment(grid::QuadGrid, A::Element, B::Element,
                       Ψs::Matrix{T}, R::T;
                       verbose=true) where {T<:Real}
    # DEBUG: only for closed shell... Otherwise you have to track the alpha and beta part..
    (verbose) && (@warn "For closed shell only!")
    ρ = 2*sum(map(x->x .^2, eachcol(Ψs))) # density associated to Ψ

    dipoles = zeros(3)
    dipoles[end] = (R/2)*(B.charge - A.charge) # nuclear contribution in the zz dipole moment
    for i in 1:3
        dipole_i = map(enumerate(grid.points)) do (ir, r)
            grid.weights[ir]*r[i]*ρ[ir]
        end
        dipoles[i] -= sum(dipole_i)
    end
    dipoles
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

function polarisability(mol::PyObject, Ez::T; FD_type=:four_point,
                        ε=1e-8, kwargs...) where {T<:Real}
    @assert FD_type ∈ (:two_point, :four_point)

    # Extract data from the mol object
    A = mol.atom_symbol(0) # element 1
    B = mol.atom_symbol(1) # element 2
    R = norm(mol.atom_coord(1) - mol.atom_coord(0))
    helfem_energy(ez) = run_helfem(A, B, R; write_hdf5=false, Ez=Ez+ez, kwargs...).e_tot
    
    # Run Helfem and computing finite difference of the energy with constant
    # electronic field in the z direction.
    αz = zero(R)
    if FD_type==:two_point
        Em = helfem_energy(Ez - ε/2)
        Ep = helfem_energy(Ez + ε/2)
        αz =  (Ep - Em)/ε
    else
        E2m = helfem_energy(Ez - 2*ε)
        Em  = helfem_energy(Ez - ε)
        Ep  = helfem_energy(Ez + ε)
        E2p = helfem_energy(Ez + 2*ε)
        αz = (-E2p + 8*Ep - 8*Em + E2m) / (12*ε)
    end
    -αz
end
