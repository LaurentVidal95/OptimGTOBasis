using OptimAOsDiatomic
using LineSearches
using Optim

# datadir = "/home/lvidal/computations/bases_optimization/HelFEM/ref_data/LiH/"
datadir="/home/vidall/computations/bases_optimization/data/test"
basis = "sto-2g"
num∫tol = 1e-7

ref_data = extract_ref_data(basis, datadir)
res = launch_optimization(ref_data; method=GradientDescent, linesearch=BackTracking(order=3));
