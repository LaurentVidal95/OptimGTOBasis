using OptimAOsDiatomic
import OptimAOsDiatomic as OAO

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")


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
