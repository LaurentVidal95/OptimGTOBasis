using JSON3

function dissociation_curve(mol::F, bases::Vector{String},
                            basis_names::Vector{String}, distances::Vector{T};
                            output_file="dissocation_curve.json", ) where
                            {T<:Real, F}
    # Compute RHF and FCI dissociation curve
    output = Dict{String, Any}()

    for (basis, basis_name) in zip(bases, basis_names)
        Es_basis = []
        for R in distances
            # Construct molecule PyObject
            mol_R = mol(basis, R)
            # Run RHF
            rhf = mol_R.RHF().run()
            @assert rhf.converged
            push!(Es_basis, rhf.e_tot)
        end
        output[basis_name] = Es_basis
    end
    output["interatomic_distances"] = distances

    # Write diss curve in JSON file
    open(io->JSON3.write(io, output, allow_inf=true), output_file, "w")
    output_file
end
