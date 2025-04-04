module OptimAOsDiatomic

using LinearAlgebra
using GenericLinearAlgebra        # apply cond or eigvals to matrices of Dual numbers
using SphericalHarmonicExpansions # Construction of AOs
using ThreadsX                    # Very basic multi-thread
using HDF5                        # Extract QuadGrid from checkfile.
using Printf                      # Cosmetic
using JuMP                        # Optimization wrapper
using Ipopt                       # Actual NL solver
using Optim

# Import pyscf globaly
using PyCall                      # Import basis_set_exchange
const pyscf = PyNULL()
function __init__()
    copy!(pyscf, pyimport("pyscf"))
end

# numerical integration
export QuadGrid
export default_spread_lim

# Spreads and contraction coeffs for each element
export Element
export extract_elements

# TODO: modify.. the AO struct is useless now
# Also create a basis struct that ease the writing and use with pyscf
# atomic orbitals
export BasisSet
include("structures/QuadGrid.jl")
include("structures/Element.jl")
include("structures/AO.jl")
include("structures/BasisSet.jl")

# Construct and solve optimization problem
export ProjectionCriterion
export EnergyCriterion
export j_proj_diatomic         # optimize ||Ψ_ref - proj(Ψ_ref)||_A
export j_E_diatomic            # optimize ||E(X) - E_ref||^2
export launch_Optim
export launch_Ipopt
include("optimize/optimization_criteria.jl")
include("optimize/setup_optim.jl")
include("optimize/optimize.jl")

# HelFEM and PySCF
export generate_basis_file
export read_helfem_data
export pyscf_property
include("external/helfem.jl")
include("external/pyscf.jl")
include("external/ipopt.jl")

export compute_pyscf_properties
export eval_criteria
export generate_molden
export plot_quantity
include("postprocess/quantities.jl")
include("postprocess/ref_quantities.jl")
include("postprocess/plot.jl")

end # module
