using JSON3

import Base.vec

struct BasisSet
    name
    str
    tag
end
BasisSet(name::String) = BasisSet(name, name, :standard)
tag(B::BasisSet) = B.tag

function Base.vec(B::BasisSet)
    split_str = split.(split(B.str, "\n"; keepempty=false), " "; keepempty=false)
    el = split_str[1][1]

    # remove label of orbitals to extract only numbers
    numbers = map(str->parse.(Float64, str), filter(x->el∉x, split_str))

    # Parse exps
    exps = [x[1] for x in numbers]
    # Parse coeffs by shells
    multipop!(tab, N) = [popfirst!(tab) for x in tab[1:N]]
    i_shell = findall(x->el ∈ x, split_str)
    len_shells = [i_shell[j+1] - (i_shell[j]+1) for j in 1:length(i_shell)-1]
    append!(len_shells, length(split_str) - i_shell[end])

    coeffs = []
    for len in len_shells
        shell = multipop!(numbers, len)
        for i_ctr in 2:length(shell[1])
            append!(coeffs, [x[i_ctr] for x in shell])
        end
    end

    Vector{Float64}(vcat(exps, coeffs))
end

function read_basis_file(Z_el::Int, ref_basis::String, file::String)
    @assert isfile(file)
    data = open(JSON3.read, file)
    basis_sets = BasisSet[]

    El = only(parse_bse_elements([Z_el], ref_basis))
    for (basis_name, coeffs) in data
        basis_str = basis_name==:standard ? ref_basis : basis_string(Vector(coeffs), El)
        push!(basis_sets, BasisSet(ref_basis, basis_str, basis_name))
    end
    basis_sets
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
