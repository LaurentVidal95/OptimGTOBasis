function compute_pyscf_properties(basis_sets::Vector{Tuple{BasisSet, BasisSet}},
                                  interatomic_distances::Vector{T}) where {T<:Real}
    # Initialize output data
    data_out = (; interatomic_distances)

    # Compute pyscf properties on each given interatomic_distances
    local_properties = (:energy, :dipole, :quadrupole, :forces, :energy_second_derivative)

    data_basis = []
    for (basis_A, basis_B) in basis_sets
        @assert basis_A.tag == basis_B.tag
        # First compute local properties for each interatomic distance
        basis_loc_prop = [[] for i in 1:5]
        for R in interatomic_distances
            # Compute RHF ground state for current R
            mol_R = mol(basis_A, basis_B, R)
            rhf = mol_R.RHF().run()
            @assert rhf.converged

            # Compute locat properties for given R
            for (i, property) in enumerate(local_properties)
                args_prop = i<4 ? (mol_R, rhf) : (basis_A, basis_B, R)
                push!(basis_loc_prop[i], pyscf_property(property, args_prop...))
            end
        end

        # Compute non local properties
        iR_mid = div(length(interatomic_distances), 2)
        R_guess = interatomic_distances[iR_mid]
        R_eq = pyscf_property(:eq_interatomic_distance, basis_A, basis_B, R_guess)

        # Add data to the full list
        push!(data_basis, merge(NamedTuple{local_properties}(basis_loc_prop), # local properties
                                (; eq_interatomic_distance=R_eq)               # non local
                                )
              )
    end
    tags = map(Bs -> Bs[1].tag, basis_sets)
    data_out = merge(data_out, NamedTuple{Tuple(tags)}(data_basis))
    return merge(data_out, (; basis=basis_sets[1][1].name))
end

function eval_criteria(basis_sets::Vector{Tuple{BasisSet, BasisSet}}, ref_data)
    A, B = ref_data.elements
    Ecrit = EnergyCriterion(ref_data)
    P0crit = ProjectionCriterion(ref_data)
    P1crit = ProjectionCriterion(ref_data; norm_type=:H¹)

    data = []
    n_data = length(ref_data.interatomic_distances)
    for (basis_A, basis_B) in basis_sets
        @assert basis_A.tag == basis_B.tag
        @info "Computing criteria for basis $(basis_A.tag)"
        data_E  = map(i -> Ecrit(basis_A, basis_B,  i),  1:n_data)
        data_L2 = map(i -> P0crit(basis_A, basis_B, i), 1:n_data)
        data_H1 = map(i -> P1crit(basis_A, basis_B, i), 1:n_data)
        push!(data, (; energy=data_E, L2_projection=data_L2,
                      H1_projection=data_H1))
    end
    tags = map(Bs -> Bs[1].tag, basis_sets)
    NamedTuple{Tuple(tags)}(data)
end
    
