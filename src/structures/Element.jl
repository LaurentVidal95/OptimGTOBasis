import Base.vec

# Just renaming
shelltype(T) = NamedTuple{(:exps, :coeffs), Tuple{Vector{T}, Matrix{T}}}

"""
Simple struct to organize all parameters by shells
"""
struct Element{T<:Real}
    name         ::String
    charge       ::Int
    basis        ::String
    shells       ::Vector{shelltype(T)}
    # Shape of the exps vectors and coeffs matrices to switch from
    # single vector to Element form.
    shape_exps   ::Vector{Int64}
    shape_coeffs ::Vector{Tuple{Int64, Int64}}
end
function Element(name::String, charge::Int, basis::String,
                 shells::Vector{shelltype(T)}) where {T<:Real}
    shape_exps = [length(shell.exps) for shell in shells]
    shape_coeffs = [size(shell.coeffs) for shell in shells]
    Element(name, charge, basis, shells, shape_exps, shape_coeffs)
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
function Element(X_vec::Vector{T1}, X_ref::Element{T2};
                 exponentiate_exps=false) where {T1, T2 <:Real}
    N_exps = len_exps(X_ref)
    exps_vec = exponentiate_exps ? exp.(X_vec[1:N_exps]) : X_vec[1:N_exps]
    coeffs_vec = X_vec[N_exps+1:end]

    pop_many!(x,n) = [popfirst!(x) for _ in 1:n]

    shells = shelltype(T1)[]
    for (α, β) in zip(X_ref.shape_exps, X_ref.shape_coeffs)
        exps = pop_many!(exps_vec, α)
        coeffs = reshape(pop_many!(coeffs_vec, prod(β)), β)
        push!(shells, (; exps, coeffs))
    end
    Element(X_ref.name, X_ref.charge, X_ref.basis, shells)
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
    Z1, Z2 = Int64(data["Z1"]), Int64(data["Z2"])
    # Extract basis data for given element using basis_set_exchange
    bse_elements = parse_bse_elements([Z1, Z2], basis)

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
function parse_bse_elements(element_charges::Vector{Int}, basis_name::String)
    element_names = element_name.(element_charges)
    bse = pyimport("basis_set_exchange")
    basis = bse.get_basis(basis_name, element_names, make_general=true)
    parsed_shells = parse_bse_shells.(values(basis["elements"]))
    # Return everything in Element structure
    [Element(name, charge, basis_name, shells) for (name, charge, shells) in
                               zip(element_names, element_charges, parsed_shells)]
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
