import Base.vec
import Base.show

# Just renaming
shelltype(T) = NamedTuple{(:exps, :coeffs), Tuple{Vector{T}, Matrix{T}}}

"""
Simple struct to organize all parameters by shells
"""
struct Element{T<:Real}
    name::String
    basis::String
    shells::Vector{shelltype(T)}
    # Shape of the exps vectors and coeffs matrices to switch from
    # single vector to Element form.
    shape_exps::Vector{Int64}
    shape_coeffs::Vector{Tuple{Int64, Int64}}
end
function Element(name::String, basis::String, shells::Vector{shelltype(T)}) where {T<:Real}
    shape_exps = [length(shell.exps) for shell in shells]
    shape_coeffs = [size(shell.coeffs) for shell in shells]
    Element(name, basis, shells, shape_exps, shape_coeffs)
end

"""
The two following function allow to pass to the storage of
exponents and coeffs in shells as storage in a single vector.
"""
# shells -> vector
function vec(X::Element{T}) where {T<:Real}
    exps = vcat([shell.exps for shell in X.shells]...)
    coeffs = vcat([vec(shell.coeffs) for shell in X.shells]...)
    vcat(exps, coeffs)
end

len_exps(X::Element{T}) where {T<:Real} = length(vcat([shell.exps for shell in X.shells]...))

# vector -> shells
# The ForwardDiff type T1 has to be different from T2
function Element(X_vec::Vector{T1}, X_ref::Element{T2}) where {T1, T2 <:Real}
    N_exps = sum(X_ref.shape_exps)
    exps_vec = X_vec[1:N_exps]
    coeffs_vec = X_vec[N_exps+1:end]

    pop_many!(x,n) = [popfirst!(x) for _ in 1:n]

    shells = shelltype(T1)[]
    for (α, β) in zip(X_ref.shape_exps, X_ref.shape_coeffs)
        exps = pop_many!(exps_vec, α)
        coeffs = reshape(pop_many!(coeffs_vec, prod(β)), β)
        push!(shells, (; exps, coeffs))
    end
    Element(X_ref.name, X_ref.basis, shells)
end

"""
Construct the Element objects associated to a given HelFEM file and GTO basis.
Output:
   • Two vectors "A" and "B"  containing the respective contracting coeffs and
     exponent of two elements for the given basis.
     Each vector is composed of NamedTuples with args "exps" and "coeffs".
     Each NamedTuple corresponds to a shell (s, p, d, ...)
"""
function extract_elements(data::Dict{String, Any}, basis::String)
    # Extract info from data
    Z1, Z2 = data["Z1"], data["Z2"]
    elements = element_name.([Z1,Z2])
    # Create GTO basis coefficient and exponent for each elements
    # using basis_set_exchange
    bse_elements = parse_bse_elements(elements, basis)

    # Case Z1=Z2
    if length(bse_elements) == 1
        Z = only(bse_elements)
        return Z, Z
    end
    return bse_elements
end

function extract_elements(datafile::String, basis::String)
    file = h5open(datafile)
    data = read(file)
    close(file)
    extract_elements(data, basis)
end

"""
Gives the name of an element given its position in the periodic table
(this data is given by HelFEM)
"""
function element_name(Z::T) where {T}
    @assert(Z<19)
    names = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
             "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
    periodic_table = Dict(T(i)=>names[i] for i in 1:18)
    periodic_table[Z]
end

"""
Return a vector whose elements contain a list of shell
and associated exponent and coefficients of contracted Gaussians
for that shell.
"""
function parse_bse_elements(elements_names::Vector{String}, basis_name::String)
    bse = pyimport("basis_set_exchange")
    basis = bse.get_basis(basis_name, elements_names, make_general=true)
    parsed_shells = parse_bse_shells.(values(basis["elements"]))
    [Element(name, basis_name, shells) for (name, shells) in zip(elements_names, parsed_shells)]
end

"""
Return exponent and corresponding contracting coefficients of a given element
(as a PyObject given by basis_set_exchange) for each shell (s, p, ...)
"""
function parse_bse_shells(element; T=Float64)
    shells = element["electron_shells"]
    parsed_shells = shelltype(T)[]
    for shell in shells
        exps = parse.(T, shell["exponents"])
        coeffs = parse.(T, shell["coefficients"])'
        push!(parsed_shells, (;exps, coeffs))
    end
    parsed_shells
end

"""
From given elements and elements name, write an AO basis file in NWChem format
(the one that seems closer to our data structure and that is understood by pyscf).
"""
function basis_string(Elements::Vector{Element{T}}) where {T<:Real}
    # Loop over all Elements
    basis_str = ""
    for El in Elements
        basis_str *= basis_string(El)
    end
    basis_str
end
function basis_string(El::Element{T}) where {T<:Real}
    basis_str = ""
    shells_names = ["S","P","D","F","G", "H"]

    for (i, shell) in enumerate(El.shells)
        shell_name = shells_names[i]
        basis_str *= "$(El.name)   $(shell_name)\n"
        # header
        mat2write = hcat(shell.exps, shell.coeffs)
        n_AO = size(shell.coeffs,2)
        for row in eachrow(mat2write)
            fmt =  Printf.Format("     "*"%10.8f    "^(n_AO+1))
            basis_str *= Printf.format(fmt, row...)*"\n"
        end
    end
    basis_str
end
basis_string(X::Vector{T}, El::Element{T}) where {T<:Real} =
    basis_string(Element(X, El))

function save_basis(Elements::Vector{Element{T}}, file) where {T<:Real}
    basis_str = basis_string(Elements)
    open(file, "w") do output_file
        println(output_file, basis_str)
    end
    nothing
end

function Base.show(io::IO, X::Element)
    println(io, "Element: $(X.name)")
    println(io, "Basis type: $(X.basis)")
    println(io, basis_string([X]))
end

function Element(El_name::String, basis::String, coeffs::Vector{T}) where {T<:Real}
    X_ref = only(extract_coeffs_and_exponents([El_name], basis))
    Element(coeffs, X_ref)
end
