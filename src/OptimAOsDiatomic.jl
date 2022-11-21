module OptimAOsDiatomic

using PyCall                      # import basis_set_exchange
using HDF5                        # Extracting QuadGrid from checkfile.
using SphericalHarmonicExpansions # Construction of AOs
using Printf                      # Printing AO labels
using LinearAlgebra
using ThreadsX

export QuadGrid
export default_spread_lim
export Element
export generate_basis_file
export construct_AOs
export eval_AOs
export extract_elements
export reference_eigenvectors
export j_L2_diatomic
export setup_bounds
include("numerical_integration.jl")
include("Element_struct.jl")
include("AO_struct.jl")
include("projection_criterion.jl")

end # module
