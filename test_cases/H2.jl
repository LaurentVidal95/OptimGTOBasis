using OptimAOsDiatomic
using JuMP

# Choose solver
using Ipopt

datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-3]...,"data/H2")
basis = "cc-pvdz"
num∫tol = 1e-9

# H2_data = extract_ref_data(basis, datadir)
H2_data_eq = extract_ref_data(basis, [joinpath(datadir,"H2_eq.hdf5")])
model = setup_optim_model(H2_data; num∫tol)
# set_optimizer_attribute(model, "max_iter", 2)
optimize!(model)
