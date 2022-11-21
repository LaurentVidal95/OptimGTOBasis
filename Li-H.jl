using OptimAOsDiatomic
using Optim
using HDF5

# extract data file
file = h5open("data/H-F/data_1.hdf5")
data = read(file)
close(file)

# Construct elements and grid
basis = "sto-3G"
H, RH, Li, RLi =  extract_elements(data, basis)
grid = QuadGrid(data)

# Extract eigenvectors and check normalization
ΨA, ΨB = reference_eigenvectors(data)

# # Launch optimization
X_init = vcat(vec(H), vec(Li))
nH = length(vec(H))
int_tol=1e-7
upper, lower = setup_bounds(grid, H, Li; int_tol)

function launch_optimization()
    @info "GTO basis optimization\n"*
        "basis: $(basis)\n"* "ζ_max: $(upper[1])\n"*
        "Maximum numerical integration error: $(int_tol)"
    if maximum(X_init) > upper[1]
        error("Some GTO cannot be integrated on the grid with"*
              " an error inferior to $(int_tol). Lower int_tol at your own risks.")
    end
    res = optimize(X->j_L2_diatomic(X, H, Li, RH, RLi, ΨA, grid),
                   lower, upper, X_init, Fminbox(LBFGS()),
                   Optim.Options(show_trace=true, extended_trace=true), autodiff=:forward);
end
