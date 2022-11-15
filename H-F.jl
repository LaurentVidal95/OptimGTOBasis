using OptimAOsDiatomic
using Optim
using HDF5

# extract data file
file = h5open("data/H-F/data_1.hdf5")
data = read(file)
close(file)

# Construct elements and grid
basis = "sto-2G"
H, RH, F, RF =  extract_elements(data, basis)
grid = QuadGrid(data)

# Extract eigenvectors and check normalization
ΨA, ΨB = reference_eigenvectors(data)

# Launch optimization
X_init = map(x-> iszero(x) ? 1e-5 : x, vcat(vec(H), vec(F)))
nH = length(vec(H))

maximum_ζ(A::Element) = maximum(maximum.([shell.exps for shell in A.shells]))
maximum_c(A::Element) = maximum(maximum.([shell.coeffs for shell in A.shells]))

"""
The optimization has to be constrained to positive parameters. As a result
the X in entry is the log of the actual parameters to so that exp(X) is always positive.
"""
function j2opt(X_log::Vector{T};
               spread_lim=default_spread_lim(grid; int_tol=1e-13),
               # default value for ctr_coeff_lim associated to conditioning to be found
               ctr_coeff_lim=1e2) where {T<:Real}
    # Retrieve positive parameters
    X = exp.(X_log)

    # Reshape the vectors XA and XB as shells to be understood by construct_AOs
    XA, XB = X[1:nH], X[nH+1:end]
    A = Element(XA, H)
    B = Element(XB, F)
    # Checks that the primitive with minimum spread can be itegrated
    # and that the contracting coefficient is not to high
    ζ_max = maximum([maximum_ζ(A), maximum_ζ(B)])
    c_max = maximum([maximum_c(A), maximum_c(B)])
    ((ζ_max > spread_lim) | (c_max > ctr_coeff_lim)) && (return Inf)
    
    # return j to minimize
    j_diatomic(A, B, RH, RF, ΨA, grid)
end


# function prune_zeros(X)
#     ids = [x[1] for x in filter(x->iszero(x[2]), collect(enumerate(X)))]
#     non_zero_ids = filter(x->x∉ ids, 1:length(X))
#     X[non_zero_ids], ids
# end
# X_init, ids = prune_zeros(vcat(vec(H), vec(F)))

# Uncomment and use X_init with prune_zeros to only optimize non zero coeffs
# N = length(X_log) + length(ids)
# X = zeros(T, N)
# X[filter(x->x∉ ids, 1:N)] .= exp.(X_log)
