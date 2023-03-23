using OptimAOsDiatomic
using JuMP

# Choose solver
using Ipopt

datadir="/home/vidall/computations/bases_optimization/data/Li-H/three_samples"
basis = "cc-pvdz"
num∫tol = 1e-10

ref_data = extract_ref_data(basis, datadir)
# model = setup_optim_model(ref_data, Ipopt; num∫tol)
# set_optimizer_attribute(model, "max_iter", 2)
# optimize!(model)
