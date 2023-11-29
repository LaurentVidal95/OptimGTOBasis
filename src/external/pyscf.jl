"""
From given elements and elements name, write an AO basis file in NWChem format
(the one that seems closer to our data structure and that is understood by pyscf).
"""
function generate_basis_file(Elements::Vector{Element{T}}, El_names::Vector{String}, file) where {T<:Real}
    @assert(length(Elements)==length(El_names))
    shells_names = ["S","P","D","F","G"]

    # Loop over all Elements
    open(file, "w") do fb
        for (El, El_name) in zip(Elements, El_names)
            for (i, shell) in enumerate(El.shells)
                shell_name = shells_names[i]
                println(fb, "$(El_name)   $(shell_name)")
                # header
                mat2write = hcat(shell.exps, shell.coeffs)
                # for AO in eachcol(shell.coeffs)                    
                #     mat2write = hcat(shell.exps, AO)
                n_AO = size(shell.coeffs,2)
                for row in eachrow(mat2write)
                    fmt =  Printf.Format("     "*"%10.8f    "^(n_AO+1))
                    println(fb, Printf.format(fmt, row...))
                end
                # end
            end
        end
    end
    nothing
end

function dissociation_curve(mol::F, basis::String, basis_opti::String, distances::Vector{T};
                            samples=[],
                            output="dissocation_curve.json",
                            FCI=false) where {T<:Real, F}
    # Compute RHF and FCI dissociation curve
    E_std, E_opti, E_fci_std, E_fci_opti = [], [], [], []

    for R in distances
        # Construct molecule PyObject
        mol_std = mol(basis, R)
        mol_opti = mol(basis_opti, R)
        # Run RHF
        rhf = mol_std.RHF().run()
        rhf_opti = mol_opti.RHF().run()
        push!(E_std, rhf.e_tot)
        push!(E_opti, rhf_opti.e_tot)
        # Run FCI
        if FCI
            fci = pyscf.fci.FCI(rhf)
            fci_opti = pyscf.fci.FCI(rhf_opti)
            push!(E_fci_std, fci.kernel()[1])
            push!(E_fci_opti, fci_opti.kernel()[1])
        end
    end

    # Write diss curve in JSON file
    output["RHF_std"] = E_std
    output["RHF_opt"] = E_opt
    output["FCI_std"] = E_std
    output["FCI_opt"] = E_opt
    output["basis"] = basis
    output["Samples"] = samples

    open(io->JSON3.write(io, output, allow_inf=true), output, "w")
    nothing
end
