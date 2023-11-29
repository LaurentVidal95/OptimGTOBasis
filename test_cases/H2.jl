using OptimAOsDiatomic
using JuMP

# Optimization
using Ipopt

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")
basis = "cc-pvdz"
num∫tol = 1e-9

function optimize_AO_basis(basis::String, datadir::String;
                           maxiter=50)
    H2_data = extract_ref_data(basis, datadir)#[joinpath(datadir,"H2_eq.hdf5")])
    model = setup_optim_model(H2_data; num∫tol)
    set_optimizer_attribute(model, "max_iter", maxiter)
    optimize!(model)
    model
end

# Computation of the dissociation curve
H2_mol(basis::String, R; symmetry=false) =
    pyscf.M(
        atom = "H 0.0 0.0 -$R/2;
                H 0.0 0.0 $R/2",
        basis=basis, # A modifier symmetry,
        symmetry=symmetry,
        unit="bohr",
        spin=0,
        charge=0)

LiH_mol(basis::String, R; symmetry=false) =
    pyscf.M(
        atom = "H 0.0 0.0 -$R/2;
                Li 0.0 0.0 $R/2",
        basis=basis, # A modifier symmetry,
        symmetry=symmetry,
        unit="bohr",
        spin=0,
        charge=0)
