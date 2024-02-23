using Plots

function plot_quantity(ref_data, pyscf_data, quantity;
                       R_eq_ref=nothing,
                       plot_error = false,
                       plot_kwargs...)
    tags = [:standard, :energy, :L2_projection, :H1_projection]
    p = plot()
    if !(plot_error)
        plot!(p, ref_data.interatomic_distances, ref_data[quantity];
              label="helFEM", linecolor=:black)
    end
    
    for tag in tags
        plot_data = pyscf_data[tag][quantity]
        !isone(length(plot_data[1])) && (plot_data = map(x->x[3], plot_data))
        if plot_error
            plot_data = abs.(plot_data - ref_data[quantity])
        end
        # Select only the z component for vector valued quantity
        plot!(p, pyscf_data.interatomic_distances, plot_data;
              label="$(tag)", markershape=:cross, markerstrokealpha=0)
    end

    # Add title etc
    xlabel = "interatomic distance (a.u)"
    ylabel = "$(quantity) (a.u.)"
    title = "$(quantity) along dissociation"*
        "$(ref_data.system) in $(pyscf_data.basis) basis"
    if plot_error
        ylabel="abs($(quantity) - ref) (log10 a.u.)"
        plot!(p, yscale=:log10)
        title = "Error of "*title
    end

    plot!(p; title, xlabel, ylabel)
    plot!(p; plot_kwargs...)
    p
end

function visualize_density(io::String, basis_A::BasisSet, basis_B::BasisSet, R)
    # Generate molden with pyscf
    mol_R = mol(basis_A, basis_B, R)
    No, Nd = mol_R.nelec;
    Nb = convert(Int64, mol_R.nao);
    
    # Add zeros to the projection if needed
    rhf = mol_R.RHF().run()
    ρ = rhf.make_rdm1()
    pyscf.tools.cubegen.density(mol_R, io, ρ)
    nothing
end

function visualize_density(io::String, ref_data)
    nothing
end
