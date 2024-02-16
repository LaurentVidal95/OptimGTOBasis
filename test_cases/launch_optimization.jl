using OptimAOsDiatomic
import OptimAOsDiatomic as OAO

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")

function optimize_AO_basis(basis::String, datadir::String, criterion_type;
                           optimizer=:Optim, # Choose between Optim.jl and JuMP (Ipopt)
                           maxiter=50,
                           norm_type=:L²,
                           guess=:bse,
                           kwargs...)
    ref_data = extract_ref_data(basis, datadir)
    criterion = criterion_type(ref_data; norm_type)
    X_guess = OAO.set_starting_point(ref_data; guess)

    if (optimizer==:Ipopt)
        return launch_Ipopt(ref_data, criterion, X_guess; maxiter, kwargs...)
    else
        return launch_Optim(ref_data, criterion, X_guess; maxiter, kwargs...)
    end
    error("Not supposed to happen")
end

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
                
