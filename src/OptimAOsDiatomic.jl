module OptimAOsDiatomic

using PyCall                      # import basis_set_exchange
using HDF5                        # Extracting QuadGrid from checkfile.
using SphericalHarmonicExpansions # Construction of AOs
using Printf                      # Printing AO labels
using LinearAlgebra
using ThreadsX
using Optim

export QuadGrid
export default_spread_lim
include("numerical_integration.jl")

export Element
export extract_elements
export generate_basis_file
include("Element_struct.jl")

export construct_AOs
export eval_AOs
include("AO_struct.jl")

export j_L2_diatomic
export reference_eigenvectors
export setup_bounds
export extract_ref_data
export launch_optimization
include("projection_criterion.jl")
include("setup_optim.jl")

end # module
