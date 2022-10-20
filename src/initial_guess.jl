using PyCall

"""
Return a vector of elements. Each element contains a list of shell
and associated exponent and coefficients of contracted Gaussians for that shell.
"""
function initial_guess(elements::Vector{String}, basis_name::String)
    bse = pyimport("basis_set_exchange")
    basis = bse.get_basis(basis_name, elements, make_general=true)
    parse_element.(values(basis["elements"]))
end

"""
Return exponent and corresponding contracting coefficients of a given element
(as a PyObject given by basis_set_exchange) for each shell (s, p, ...)
"""
function parse_element(element; T=Float64)
    shells = element["electron_shells"]
    parsed_shells = NamedTuple{(:exps, :coeffs),
                               Tuple{Vector{T}, Matrix{T}}}[]
    for shell in shells
        exps = parse.(T, shell["exponents"])
        coeffs= parse.(T, shell["coefficients"])'
        push!(parsed_shells, (;exps, coeffs))
    end
    parsed_shells
end
