using OptimAOsDiatomic
import OptimAOsDiatomic as OAO
using JuMP
using PyCall

# Optimization
using Ipopt

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")
num∫tol = 1e-9

function optimize_AO_basis(basis::String, datadir::String;
                           maxiter=50)
    H2_data = extract_ref_data(basis, datadir)#[joinpath(datadir,"helfem_1.4.hdf5")])
    model = setup_optim_model(H2_data; num∫tol)
    set_optimizer_attribute(model, "max_iter", maxiter)
    optimize!(model)
    model
end

# Computation of the dissociation curve
# pyscf = pyimport("pyscf")

H2_mol(basis::String, R; symmetry=false) =
    OAO.pyscf.M(;
        atom = "H 0.0 0.0 0.0;
                H 0.0 0.0 $R",
        basis,
        symmetry,
        unit="bohr",
    )
LiH_mol(basis::String, R; symmetry=false) =
    OAO.pyscf.M(;
        atom = "H 0.0 0.0 0.0;
                Li 0.0 0.0 $R",
        basis,
        symmetry,
        unit="bohr",
    )
