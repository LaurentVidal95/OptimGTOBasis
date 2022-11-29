"""
From given elements and elements name, write an AO basis file in NWChem format
(the one that seems closer to our data structure and that is understood by pyscf).
"""
function generate_basis_file(Elements::Vector{Element{T}}, El_names::Vector{String}, file) where {T<:Real}
    @assert(length(Elements)==length(El_names))
    shells_names = ["S","P","D","F","G"]

    # Loop over all Elements
    open(file, "w") do fb
        for (El, El_name) in zip(Elements, El_names)
            for (i, shell) in enumerate(El.shells)
                shell_name = shells_names[i]
                # header
                for AO in eachcol(shell.coeffs)
                    println(fb, "$(El_name)   $(shell_name)")
                    mat2write = hcat(shell.exps, AO)
                    for row in eachrow(mat2write)
                        @printf fb "     %-10.8f %-10.8f \n" row...
                    end
                end
            end
        end
    end
    nothing         
end

# TODO: Launch pyscf RHF computation
