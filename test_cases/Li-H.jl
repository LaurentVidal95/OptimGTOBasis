using OptimAOsDiatomic
using JuMP

datadir="/home/vidall/computations/bases_optimization/data/test/"
basis = "sto-2g"
num∫tol = 1e-7

ref_data = extract_ref_data(basis, datadir)
# model = setup_optim_model(ref_data; num∫tol=1e-8)
# optimize!(model)
