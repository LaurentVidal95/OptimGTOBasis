module OptimAOsDiatomic

using LinearAlgebra
using GenericLinearAlgebra        # apply cond or eigvals to matrices of Dual numbers
using SphericalHarmonicExpansions # Construction of AOs
using ThreadsX                    # Very basic multi-thread
using HDF5                        # Extract QuadGrid from checkfile.
using Printf                      # Cosmetic
using JuMP                        # Optimization wrapper
using Ipopt                       # Actual NL solver

# Import pyscf globaly
using PyCall                      # Import basis_set_exchange
const pyscf = PyNULL()
function __init__()
    copy!(pyscf, pyimport("pyscf"))
end

# # DEBUG
# using LegendrePolynomials
# using Optim
#

# numerical integration
export QuadGrid
export default_spread_lim

# Spreads and contraction coeffs for each element
export Element
export extract_elements

# atomic orbitals
export construct_AOs
export eval_AOs
include("structures/QuadGrid.jl")
include("structures/Element.jl")
include("structures/AO.jl")

# Extract data from basis_set_exchange
export reference_eigenvectors
export extract_ref_data

# Construct and solve optimization problem
export j_L2_diatomic           # optimize ||Ψ_ref - proj(Ψ_ref)||_L²
export j_E_diatomic            # optimize ||E(X) - E_ref||^2
export setup_bounds!
export setup_optim_model
include("optimize/optimizatin_criteria.jl")
include("optimize/setup_optim.jl")

# HelFEM and PySCF
export generate_basis_file
include("external/helfem.jl")
include("external/pyscf.jl")
# Small useful routines
include("utils.jl")

end # module
