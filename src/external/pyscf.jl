energy(mol::PyObject, rhf::PyObject) = rhf.e_tot
function energy(mol::PyObject)
    # Check that the system is closed-shell
    No, Nd = mol.nelec; Ns = No - Nd;
    @assert iszero(Ns) "System has to be a closed shell system"
    rhf = mol.RHF().run()
    @assert rhf.converged "SCF not converged"
    energy(mol, rhf)
end

function dipole_moment(mol::PyObject, C::AbstractArray; verbose=false)
    # TODO: sanity checks
    R1 = mol.atom_coords()[1,3];  R2 = mol.atom_coords()[2,3]
    Z1 = mol.atom_charge(0);      Z2 = mol.atom_charge(1)
    μz_mat = mol.intor("int1e_z")
    ρ = 2*C*C'

    μ_elec = -tr(μz_mat*ρ)
    μ_nuc = (Z1*R1 + Z2*R2)
    if verbose
        @show "Nuclear $(μ_nuc)"
        @show "Electronic $(μ_elec)"
    end
    μ_elec + μ_nuc
end
function dipole_moment(mol::PyObject; verbose=false)
    No, Nd = mol.nelec; Ns = No - Nd;
    rhf = mol.RHF().run()
    @assert rhf.converged "SCF not converged"
    dipole_moment(mol, rhf.mo_coeff[:,No]; verbose)
end

function quadrupole_moment(mol::PyObject, C::AbstractArray; verbose=false)
    # Sanity checks
    No, Nd = mol.nelec; Ns = No - Nd;
    @assert iszero(Ns) "System has to be a closed shell system"
    @assert iszero(mol.atom_coords()[1:2,1:2]) "Atoms have to be on the Z axis"

    ρ = 2*C*C'

    Q_mat_zz = 3*mol.intor("int1e_zz") - mol.intor("int1e_r2")
    Q_mat_xx = (3/2)*(mol.intor("int1e_r2") -mol.intor("int1e_zz")) - mol.intor("int1e_r2")

    Q_elec_zz = -(1/2)*tr(Q_mat_zz*ρ)
    Q_elec_xx = -(1/2)*tr(Q_mat_xx*ρ)
    Q_elec = [Q_elec_xx, Q_elec_xx, Q_elec_zz]

    # Nuclear contribution
    Z1, Z2 = mol.atom_charge(0), mol.atom_charge(1)
    R = norm(mol.atom_coords()[1,:] - mol.atom_coords()[2,:])
    Q_nuc = [0.0, 0.0, (Z1+Z2)*(R^2)/4]

    if verbose
        @info "Nuclear: $(Q_nuc)"
        @info "Electronic: $(Q_elec)"
        @info "Total: $(Q_nuc + Q_elec)"
    end
    Q_elec + Q_nuc
end
function quadrupole_moment(mol::PyObject; verbose=false)
    No, Nd = mol.nelec; Ns = No - Nd;
    rhf = mol.RHF().run()
    @assert rhf.converged "SCF not converged"
    quadrupole_moment(mol, rhf.mo_coeff[:,No]; verbose)
end

function eq_interatomic_distance(mol::PyObject, rhf::PyObject)
    # Define python function that wrapps geometry optimization
    py"""
    def geo_opt(rhf):
        from pyscf.geomopt.geometric_solver import optimize
        mol_eq = optimize(rhf, maxsteps=100)
        return mol_eq.atom_coords()
    """
    # Launch function and extract result
    eq_coords = py"geo_opt"(rhf)
    norm(eq_coords[1,:] - eq_coords[2,:])
end
function eq_interatomic_distance(mol::PyObject)
    rhf = mol.RHF().run()
    @assert rhf.converged "SCF not converged"
    eq_interatomic_distance(mol, rhf)
end


function dissociation_curve(mol::F, bases::Vector{String},
                            basis_names::Vector{String}, distances::Vector{T};
                            output_file=nothing, ) where
                            {T<:Real, F}
    # Compute RHF and FCI dissociation curve
    output = Dict{String, Any}()

    for (basis, basis_name) in zip(bases, basis_names)
        Es_basis = []
        for R in distances
            # Construct molecule PyObject
            mol_R = mol(basis, R)
            push!(Es_basis, energy(mol_R))
        end
        output[basis_name] = Es_basis
    end
    output["interatomic_distances"] = distances

    # Write diss curve in JSON file
    !isnothing(output_file) && (open(io->JSON3.write(io, output, allow_inf=true), output_file, "w"))
    output_file
end
