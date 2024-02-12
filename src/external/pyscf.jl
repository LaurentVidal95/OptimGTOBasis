using JSON3

function dissociation_curve(mol::F, basis::String, basis_opt::String, distances::Vector{T};
                            samples=[],
                            output_file="dissocation_curve.json",
                            FCI=false) where {T<:Real, F}
    # Compute RHF and FCI dissociation curve
    E_std, E_opt, E_fci_std, E_fci_opt = [], [], [], []
    output = Dict{String, Any}()

    for R in distances
        # Construct molecule PyObject
        mol_std = mol(basis, R)
        mol_opt = mol(basis_opt, R)
        # Run RHF
        rhf = mol_std.RHF().run()
        rhf_opt = mol_opt.RHF().run()
        push!(E_std, rhf.e_tot)
        push!(E_opt, rhf_opt.e_tot)
        # Run FCI
        if FCI
            fci = pyscf.fci.FCI(rhf)
            fci_opt = pyscf.fci.FCI(rhf_opt)
            push!(E_fci_std, fci.kernel()[1])
            push!(E_fci_opt, fci_opt.kernel()[1])
        end
    end

    # Write diss curve in JSON file
    output["RHF_std"] = E_std
    output["RHF_opt"] = E_opt
    if FCI
        output["FCI_std"] = E_std
        output["FCI_opt"] = E_opt
    end
    output["bond lengths"] = distances
    output["basis"] = basis
    output["basis_opt"] = basis_opt
    output["Samples"] = samples
    open(io->JSON3.write(io, output, allow_inf=true), output_file, "w")
    output_file
end
