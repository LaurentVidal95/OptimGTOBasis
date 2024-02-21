function mol(basis_A::BasisSet, basis_B::BasisSet, R::T;
             verbose=0) where {T<:Real}
    A = basis_A.element
    B = basis_B.element
    basis = Dict([A.name => String(basis_A),
                  B.name => String(basis_B)])
    
    atom = "$(A.name) 0.0 0.0 -$(R/2);
            $(B.name) 0.0 0.0 $(R/2)"
    unit="bohr"
    pyscf.M(; atom, basis, unit, verbose)
end

"""
Wrapper around all functions bellow using pyscf.
"""
function pyscf_property(s::Symbol, args...; kwargs...)
    method = Dict(:energy => energy,
                  :dipole => dipole_moment,
                  :quadrupole => quadrupole_moment,
                  :eq_interatomic_distance => eq_interatomic_distance,
                  :forces => force_FD,
                  :energy_second_derivative=>∂2E_FD,
                  :polarizability=>nothing)
    return method[s](args...; kwargs...)
    error("This property is unaccessible for the moment")
end

energy(mol::PyObject, rhf::PyObject) = rhf.e_tot
function energy(mol::PyObject)
    # Check that the system is closed-shell
    No, Nd = mol.nelec; Ns = No - Nd;
    @assert iszero(Ns) "System has to be a closed shell system"
    rhf = mol.RHF().run()
    @assert rhf.converged "SCF not converged"
    energy(mol, rhf)
end
energy(basis_A::BasisSet, basis_B::BasisSet, R::T; verbose=0) where {T<:Real} =
    energy(mol(basis_A, basis_B, R; verbose))

function dipole_moment(mol::PyObject, rhf::PyObject; verbose=false)
    # Compute 1-RDM
    ρ = rhf.make_rdm1()

    # TODO: sanity checks
    R1 = mol.atom_coords()[1,3];  R2 = mol.atom_coords()[2,3]
    Z1 = mol.atom_charge(0);      Z2 = mol.atom_charge(1)
    μz_mat = mol.intor("int1e_z")

    μ_elec = -tr(μz_mat*ρ)
    μ_nuc = convert(Float64, (Z1*R1 + Z2*R2))
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
    dipole_moment(mol, rhf; verbose)
end
dipole_moment(basis_A::BasisSet, basis_B::BasisSet, R::T; verbose=0) where {T<:Real} =
    dipole_moment(mol(basis_A, basis_B, R; verbose); verbose)

function quadrupole_moment(mol::PyObject, rhf::PyObject; verbose=false)
    # Compute 1-RDM
    ρ = rhf.make_rdm1()

    # Sanity checks
    No, Nd = mol.nelec; Ns = No - Nd;
    @assert iszero(Ns) "System has to be a closed shell system"
    @assert iszero(mol.atom_coords()[1:2,1:2]) "Atoms have to be on the Z axis"

    Q_mat_zz = 3*mol.intor("int1e_zz") - mol.intor("int1e_r2")
    Q_mat_xx = (3/2)*(mol.intor("int1e_r2") -mol.intor("int1e_zz")) - mol.intor("int1e_r2")

    Q_elec_zz = -(1/2)*tr(Q_mat_zz*ρ)
    Q_elec_xx = -(1/2)*tr(Q_mat_xx*ρ)
    Q_elec = [Q_elec_xx, Q_elec_xx, Q_elec_zz]

    # Nuclear contribution
    Z1, Z2 = mol.atom_charge(0), mol.atom_charge(1)
    R = norm(mol.atom_coords()[1,:] - mol.atom_coords()[2,:])
    Q_nuc = convert.(Float64, [0.0, 0.0, (Z1+Z2)*(R^2)/4])

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
    quadrupole_moment(mol, rhf; verbose)
end
quadrupole_moment(basis_A::BasisSet, basis_B::BasisSet, R::T; verbose=0) where {T<:Real} =
    quadrupole_moment(mol(basis_A, basis_B, R); verbose)

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
function eq_interatomic_distance(basis_A::BasisSet, basis_B::BasisSet, R::T;
                                 verbose=0)
    eq_interatomic_distance(mol(basis_A, basis_B, R))
end

function force_FD(basis_A::BasisSet, basis_B::BasisSet,
                  R::T; ε=1e-5, verbose=false) where {T<:Real}
    E⁻ = energy(basis_A, basis_B, R - ε/2; verbose)
    E⁺ = energy(basis_A, basis_B, R + ε/2; verbose)
    (E⁺ - E⁻)/ε
end

function ∂2E_FD(basis_A::BasisSet, basis_B::BasisSet, R::T;
                ε=1e-5, verbose=false) where {T<:Real}
    F⁻ = force_FD(basis_A, basis_B, R - ε/2; verbose)
    F⁺ = force_FD(basis_A, basis_B, R + ε/2; verbose)
    (F⁺ - F⁻)/ε
end
