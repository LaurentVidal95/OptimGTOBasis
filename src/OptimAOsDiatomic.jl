module OptimAOsDiatomic

using LinearAlgebra
using SphericalHarmonicExpansions # Construction of AOs
using ThreadsX                    # Very basic multi-thread
using PyCall                      # Import basis_set_exchange
using HDF5                        # Extract QuadGrid from checkfile.
using Printf                      # Cosmetic
using JuMP                        # Optimization wrapper
using Ipopt                       # Actual NL solver

# DEBUG
using LegendrePolynomials
using Optim
#

## All basic structures of the code

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

## Optimization of the AO basis

# Extract data from basis_set_exchange
export reference_eigenvectors
export extract_ref_data

# Construct and solve optimization problem
export j_L2_diatomic # optimize ||Ψ_ref - proj(Ψ_ref)||_L²
export setup_bounds!
export setup_optim_model
include("optimize/projection_criterion.jl")
include("optimize/setup_optim.jl")

export generate_basis_file
include("compare_energies/pyscf_wrapper.jl")

# Small useful routines
include("utils.jl")
include("helfem.jl")

end # module
