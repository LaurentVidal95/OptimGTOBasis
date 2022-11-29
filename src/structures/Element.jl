import Base.vec

# Just renaming
shelltype(T) = NamedTuple{(:exps, :coeffs), Tuple{Vector{T}, Matrix{T}}}

"""
Simple struct to organize all parameters by shells
"""
struct Element{T<:Real}
    shells::Vector{shelltype(T)}
    # Shape of the exps vectors and coeffs matrices to switch from
    # single vector to Element form.
    shape_exps::Vector{Int64}
    shape_coeffs::Vector{Tuple{Int64, Int64}}
end
function Element(shells::Vector{shelltype(T)}) where {T<:Real}
    shape_exps = [length(shell.exps) for shell in shells]
    shape_coeffs = [size(shell.coeffs) for shell in shells]
    Element(shells, shape_exps, shape_coeffs)
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
    Element(shells)
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
    A, B = extract_coeffs_and_exponents(elements, basis)
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
function extract_coeffs_and_exponents(elements_names::Vector{String}, basis_name::String)
    bse = pyimport("basis_set_exchange")
    basis = bse.get_basis(basis_name, elements_names, make_general=true)
    parse_bse_element.(values(basis["elements"]))
end

"""
Return exponent and corresponding contracting coefficients of a given element
(as a PyObject given by basis_set_exchange) for each shell (s, p, ...)
"""
function parse_bse_element(element; T=Float64)
    shells = element["electron_shells"]
    parsed_shells = shelltype(T)[]
    for shell in shells
        exps = parse.(T, shell["exponents"])
        coeffs = parse.(T, shell["coefficients"])'
        push!(parsed_shells, (;exps, coeffs))
    end
    Element(parsed_shells)
end
