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
    !isnothing(output_file) && (open(io->JSON3.write(io, output, allow_inf=true), output_file, "w"))
    output_file
end

function quadrupole(mol::F, bases::Vector{String}, basis_names::Vector{String},
                    charges, distances::Vector{T};
                    ) where {T<:Real, F}
    @warn "Only for closed shell"

    # Compute RHF and FCI dissociation curve
    output = Dict{String, Any}()

    for (basis, basis_name) in zip(bases, basis_names)
        quads = []
        for R in distances
            mol_R = mol(basis, R; verbose=0)
            push!(quads, quadrupole(mol_R, charges, R))
        end
        output[basis] = quads
    end
    output
end
function quadrupole(mol::PyObject, charges, R)
    No = Int(mol.nelectron/2)
    # Compute optimal MOs
    rhf = mol.RHF().run()
    @assert rhf.converged
    C = rhf.mo_coeff[:,No]

    # Compute quadrupole operator from AOs
    Q_mat = 3*mol.intor("int1e_zz") - mol.intor("int1e_r2")
    Q_nuc = sum(charges)*(R^2)/4
    Q_elec = -tr(C'*Q_mat*C)
    @info "Nuclear: $(Q_nuc)"
    @info "Electronic: $(Q_elec)"
    @info "Total: $(Q_nuc + Q_elec)"
    # Compute moment
    Q_elec + Q_nuc
end
