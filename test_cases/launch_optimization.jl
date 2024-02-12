using OptimAOsDiatomic
import OptimAOsDiatomic as OAO

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")

function optimize_AO_basis(basis::String, datadir::String, criterion_type;
                           optimizer=:Optim, # Choose between Optim.jl and JuMP (Ipopt)
                           maxiter=50,
                           gridtol=1e-9,
                           guess=:bse,
                           kwargs...)
    ref_data = extract_ref_data(basis, datadir)
    criterion = criterion_type(ref_data; gridtol)
    X_guess = OAO.set_starting_point(ref_data; guess)

    if (optimizer==:Ipopt)
        return launch_Ipopt(ref_data, criterion, X_guess; maxiter, kwargs...)
    else
        return launch_Optim(ref_data, criterion, X_guess; maxiter, kwargs...)
    end
    error("Not supposed to happen")
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

H2O_mol(basis_H::String, basis_O::String) =
    OAO.pyscf.M(;
                atom = "O          0.00000        0.00000        0.11779;
                        H          0.00000        0.75545       -0.47116;
                        H          0.00000       -0.75545       -0.47116",
                basis=Dict(["O"=>"basis_O", "H"=>basis_H]),
                unit="angstrom"
                )
                
