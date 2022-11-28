using OptimAOsDiatomic

datadir="/home/lvidal/computations/bases_optimization/HelFEM/ref_data/test"
basis = "sto-2g"
num∫tol = 1e-7

ref_data = extract_ref_data(basis, datadir)
# res = launch_optimization(ref_data; method=GradientDescent, linesearch=BackTracking(order=3));
