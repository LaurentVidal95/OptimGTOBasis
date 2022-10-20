using PyCall

"""
Compare the overlap obtained by numerical integration with
the overlap obtained with pyscf, given the pyscf molecule
and the overlap for numerical integration.
"""
function test_ovlp(mol::PyObject, S_num::AbstractMatrix{T}) where {T<:Real}
    S_pyscf = mol.intor("int1e_ovlp")
    norm(S_pyscf - S_num)
end
function test_ovlp(atoms::Vector{String}, positions::Vector{Vector{T}},
                   basis_name::String, S_num::AbstractMatrix{T};
                   kwargs...) where {T<:Real}
    # Ensure that spin charge and units are specified to avoid bugs
    @assert :spin in keys(kwargs)   # see pyscf doc for def
    @assert :charge in keys(kwargs) # electric charge
    @assert :unit in keys(kwargs)   # bohr or angstrom
    
    pyscf = pyimport("pyscf")

    # Assemble pyscf molecule
    atom_pyscf = ""
    for (atom, position) in zip(atoms, positions)
        atom_pyscf *= atom*" $(position[1]) $(position[2]) $(position[3]); "
    end
    mol = pyscf.M(;atom=atom_pyscf, basis=basis_name, symmetry=true, kwargs...)
    # Compare both overlaps
    test_ovlp(mol, S_num)
end
