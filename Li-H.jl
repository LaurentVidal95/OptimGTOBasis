using OptimAOsDiatomic

# datadir = "/home/lvidal/computations/bases_optimization/HelFEM/ref_data/LiH/"
datadir="/home/lvidal/computations/bases_optimization/HelFEM/ref_data/test"
basis = "sto-3g"
num∫tol = 1e-7

ref_data = extract_ref_data(basis, datadir)
res = launch_optimization(ref_data);
