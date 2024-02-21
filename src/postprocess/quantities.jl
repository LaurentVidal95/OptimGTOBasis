function compute_pyscf_properties(basis_file_A::String, basis_file_B::String,
                                  interatomic_distances::Vector{T}) where {T<:Real}
    # Initialize output data
    data_out = (; interatomic_distances)

    # Compute pyscf properties on each given interatomic_distances
    local_properties = (:energy, :dipole, :quadrupole, :forces, :energy_second_derivative)

    data_basis = []
    for basis in basis_sets
        # First compute local properties for each interatomic distance
        basis_loc_prop = [[] for i in 1:5]
        for R in interatomic_distances
            # Compute RHF ground state for current R
            mol_R = mol(basis.str, R)
            rhf = mol_R.RHF().run()
            @assert rhf.converged

            # Compute locat properties for given R
            for (i, property) in enumerate(local_properties)
                args_prop = i<4 ? (mol_R, rhf) : (mol, basis.str, R)
                push!(basis_loc_prop[i], pyscf_property(property, args_prop...))
            end
        end

        # Compute non local properties
        iR_mid = div(length(interatomic_distances), 2)
        R_guess = interatomic_distances[iR_mid]
        R_eq = pyscf_property(:eq_interatomic_distance, mol, basis.str, R_guess)

        # Add data to the full list
        push!(data_basis, merge(NamedTuple{local_properties}(basis_loc_prop), # local properties
                                (; eq_interatomic_distance=R_eq)               # non local
                                )
              )
    end
    data_out = merge(data_out, NamedTuple{Tuple(tag.(basis_sets))}(data_basis))
    return merge(data_out, (; basis=basis_sets[1].name))
end

function eval_criteria(basis_sets::Vector{BasisSet}, ref_data)
    A, B = ref_data.elements
    Ecrit = EnergyCriterion(ref_data)
    P0crit = ProjectionCriterion(ref_data)
    P1crit = ProjectionCriterion(ref_data; norm_type=:H¹)

    data = []
    n_data = length(ref_data.interatomic_distances)
    for basis in basis_sets
        @info "Computing criteria for basis $(basis.tag)"
        (basis.tag==:standard) && (basis = BasisSet(basis.name, basis_string([A]), :standard))
        data_E  = map(i -> Ecrit(basis, A, B, i),  1:n_data)
        data_L2 = map(i -> P0crit(basis, A, B, i), 1:n_data)
        data_H1 = map(i -> P1crit(basis, A, B, i), 1:n_data)
        push!(data, (; energy=data_E, L2_projection=data_L2,
                      H1_projection=data_H1))
    end
    NamedTuple{Tuple(tag.(basis_sets))}(data)
end
    
