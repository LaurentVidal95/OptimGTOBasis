using OptimAOsDiatomic
import OptimAOsDiatomic as OAO
using JuMP
using Ipopt

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")

function optimize_AO_basis(basis::String, datadir::String, criterion_type;
                           maxiter=50,
                           gridtol=1e-9,
                           kwargs...)
    ref_data = extract_ref_data(basis, datadir)
    criterion = criterion_type(ref_data; gridtol)
    model = setup_optim_model(ref_data, criterion; kwargs...)
    set_optimizer_attribute(model, "line_search_method", "filter")
    set_optimizer_attribute(model, "max_iter", maxiter)
    optimize!(model)
    model
end

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
