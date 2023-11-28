using OptimAOsDiatomic
using JuMP

# Choose solver
using Ipopt

# datadir=joinpath(splitpath(pathof(OptimAOsDiatomic))[1:end-2]...,"test_cases/Li-H/data")
basis = "cc-pvdz"
num∫tol = 1e-10

# LiH_data = extract_ref_data(basis, datadir)
# model = setup_optim_model(LiH_data, Ipopt; num∫tol)
# set_optimizer_attribute(model, "max_iter", 2)
# optimize!(model)
