using OptimAOsDiatomic
using HDF5

# extract data file
file = h5open("data/H-FE/data_1.hdf5")
data = read(file)
close(file)

# Construct elements and grid
basis = "cc-pvdz"
H, RH, F, RF =  extract_elements(data, basis)
grid = QuadGrid(data)

# Extract eigenvectors and check normalization
ΨA, ΨB = reference_eigenvectors(data)
SA =  dot(grid, ΨA, ΨA)
SB =  dot(grid, ΨB, ΨB)
@show (norm(SA - I), norm(SB - I))

# Contruct AOs for H and eval on the grid
AOs_H = construct_AOs(H; position=RH, mmax=grid.mmax)
AOs_on_grid = eval_AOs(grid, AOs_H)

# Compute overlap
S = overlap(grid, AOs_on_grid)
