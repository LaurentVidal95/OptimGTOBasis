using OptimAOsDiatomic
import OptimAOsDiatomic as OAO

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")

function optimize_AO_basis(basis::String, datadir::String, criterion_type;
                           optimizer=:Optim, # Choose between Optim.jl and JuMP (Ipopt)
                           maxiter=100,
                           norm_type=:L²,
                           guess=:bse,
                           kwargs...)
    ref_data = read_helfem_data(basis, datadir)
    criterion = criterion_type(ref_data; norm_type)
    X_guess = OAO.set_starting_point(ref_data; guess)

    if (optimizer==:Ipopt)
        return launch_Ipopt(ref_data, criterion, X_guess; maxiter, kwargs...)
    else
        return launch_Optim(ref_data, criterion, X_guess; maxiter, kwargs...)
    end
    error("Not supposed to happen")
end

function compute_pyscf_data(mol::Function, basis_sets::Vector{BasisSet},
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
        R_eq = OAO.pyscf_property(:eq_interatomic_distance, mol, basis.str, R_guess)

        # Add data to the full list
        push!(data_basis, merge(NamedTuple{local_properties}(basis_loc_prop), # local properties
                                (; eq_interatomic_distance=R_eq)               # non local
                                )
              )
    end
    data_out = merge(data_out, NamedTuple{Tuple(OAO.tag.(basis_sets))}(data_basis))
    return merge(data_out, (; basis=basis_sets[1].name))
end

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

function plot_quantity(ref_data, pyscf_data, quantity;
                       R_eq_ref=nothing,
                       plot_kwargs...)
    p = plot(ref_data.interatomic_distances, ref_data[quantity];
             label="helFEM", linecolor=:black)
    tags = [:standard, :energy, :L2_projection, :H1_projection]

    for tag in tags
        plot_data = pyscf_data[tag][quantity]
        # Select only the z component for vector valued quantity
        !isone(length(plot_data[1])) && (plot_data = map(x->x[3], plot_data))
        plot!(p, pyscf_data.interatomic_distances, plot_data;
              label="$(tag)", markershape=:cross, markerstrokealpha=0)
    end

    # Add title etc
    xlabel = "interatomic distance (a.u)"
    ylabel = "$(quantity) (a.u.)"
    title = "$(quantity) along dissociation"*
        "$(ref_data.system) in $(pyscf_data.basis) basis"
    plot!(p; title, xlabel, ylabel)
    plot!(p; plot_kwargs...)
    p
end


# List of test systems
X2_mol(el::String, basis::String, R; verbose=0) =
        OAO.pyscf.M(;
                    atom = "$(el) 0.0 0.0 -$(R/2);
                            $(el) 0.0 0.0 $(R/2)",
                    basis,
                    unit="bohr",
                    verbose
                    )

H2_mol(basis::String, R; verbose=0) =
    OAO.pyscf.M(;
        atom = "H 0.0 0.0 -$(R/2);
                H 0.0 0.0 $(R/2)",
                basis,
                unit="bohr",
                verbose
                )

LiH_mol(basis::String, R; verbose=0)=
    OAO.pyscf.M(;
        atom = "H 0.0 0.0 -$(R/2);
                Li 0.0 0.0 $(R/2)",
                basis,
                unit="bohr",
                verbose
    )

H2O_mol(basis_H::String, basis_O::String; verbose=0) =
    OAO.pyscf.M(;
                atom = "O          0.00000        0.00000        0.11779;
                        H          0.00000        0.75545       -0.47116;
                        H          0.00000       -0.75545       -0.47116",
                basis=Dict(["O"=>basis_O, "H"=>basis_H]),
                unit="angstrom",
                verbose
                )
