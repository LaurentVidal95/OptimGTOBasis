"""
    struct RadialPart(exps, coeffs)

Radial part ``R_{nl}`` of an atomic orbital ``AO(X) = Y_{lm}(X)R_{nl}(X)``, as
a contracted gaussian ``R_{nl}(X) = ∑ᵢ cᵢ exp(-ζᵢ*|X|²)``.

# Fields
    - `exps`: contains the list of exponents ``ζᵢ`` 
    - `coeffs`: list of coefficients ``cᵢ``
"""
struct RadialPart{T<:Real}
    exps   ::Vector{T}
    coeffs ::Vector{T}
end
(Rnl::RadialPart)(X) = Rnl.coeffs'exp.(map(ζ -> -ζ*norm(X)^2, Rnl.exps))

"""
    struct AO(Rnl, Ylm, center)

Generic atomic orbital (AO) basis function ``AO(X) = Y_{lm}(X)R_{nl}(X)``.

# Fields
    - `Rnl`: the radial part ``R_{nl}`` is defined in the above struct [`RadialPart`](@ref)
    - `Ylm`: real spherical harmonic, defined with the `SphericalHarmonicExpansions.jl`
      package. Note that it contains the ``r^l`` part of the atomic orbitals, so that
      it is consistent with the Radial part definition.
    - `center`: the center of the AO.
"""
struct AO{T1, T2 <: Real, F}
    Rnl    ::RadialPart{T1}
    Ylm    ::F
    center ::Vector{T2}
end
function AO(l::TI, m::TI, exps::Vector{T1}, coeffs::Vector{T1},
            center::Vector{T2}) where {T1, T2<:Real, TI <:Int}
    Rnl = RadialPart(exps, coeffs)
    @polyvar x y z
    Ylm(r1,r2,r3) = rlylm(l, m, x, y, z)((x,y,z)=>(r1,r2,r3))
    AO(Rnl, Ylm, center)
end
(Ψ::AO)(X::Vector) = Ψ.Rnl(X .- Ψ.center)*Ψ.Ylm((X .- Ψ.center)...)

"""
Element is a vector of shells, each being a named tuple with args
exps and coeffs, corresponding to the Gaussians for that shell.
"""
function AO_basis(X::Element{T};
                  position=zeros(Float64, 3), # [0., 0., ± Rh]
                  mmax,
                  verbose=true) where {T<:Real}
    AOs = AO[]
    (verbose) && (@info "Added AOs: ")
    # Run other n, l and m ∈ {-l,...,l} and add each corresponding AO
    for (i_shell, shell) in enumerate(X.shells)
        l = i_shell-1
        for n in 1:size(shell.coeffs,2)
            for m in -min(mmax, l):min(mmax,l)
                label = ["s","p","d","f","g"][i_shell]
                (verbose) && (println(@sprintf("%-10s %-6s",
                              "Χ_{$(n+l)$l$m}", "$(n+l)$label")))
                push!(AOs, AO(l, m, shell.exps, vec(shell.coeffs[:,n]), position))
            end
        end
    end
    AOs
end

function eval_abs_ms(mol::PyObject)
    # Compute AO labels for molecule X₂
    all_labels = mol.ao_labels(fmt=false)
    N_aos = div(length(all_labels),2)
    labels = all_labels[1:N_aos] # AO label on first X

    # Compute abs of the associated quantum number m
    multipopfirst!(tab, N) = [popfirst!(tab) for x in tab[1:N]]
    abs_ms = Int[]
    while !isempty(labels)
        l = split(labels[1][3],"")[end]
        n_aos_nl=Dict("s"=>1, "p"=>3, "d"=>5, "f"=>7, "g"=>9,
                      "h"=>11)
        aos_nl = multipopfirst!(labels, n_aos_nl[l])
        if l=="s"
            append!(abs_ms, [0])
        elseif l=="p"
            append!(abs_ms, [1,1,0])
        else # l=d,f,g,...
            # Compute numbers of aos with same n and l
            half_n_ao_nl = div(length(aos_nl), 2)            
            append!(abs_ms, abs.(collect(-half_n_ao_nl:half_n_ao_nl)))
        end
    end
    vcat(abs_ms, abs_ms)
end


function eval_AOs(grid::QuadGrid, A::Element{T1}, B::Element{T1},
                  R::T2) where {T1, T2<:Real}
    basis = Dict([A.name => basis_string([A]),
                  B.name => basis_string([B])])
    tmp_mol = pyscf.M(;atom="$(A.name) 0.0 0.0 0.0;
                             $(B.name) 0.0 0.0 $R",
                      basis
                      )
    X = pyscf.dft.numint.eval_ao(tmp_mol, grid.points)
    # filter X with mmax: X[:,i]=0 if |m(X[:,i])| > mmax
    @show abs_ms = eval_abs_ms(tmp_mol)
    for (im, m) in enumerate(abs_ms)
        if m>grid.mmax
            X[:,im] .= zero(X[:,im])
        end
    end
    X
end

function overlap(X::Element{T}, R; norm_type=:L²) where {T<:Real}
    basis = Dict([X.name => basis_string([X])])
    tmp_mol = pyscf.M(;atom="$(X.name) 0.0 0.0 0.0;
                             $(X.name) 0.0 0.0 $R",
                      basis
                      )
    (norm_type=:L²) && (return tmp_mol.intor("int1e_ovlp"))
    2*tmp_mol.intor("int1e_kin")
end

eval_AO(grid::QuadGrid, AO::AO) = ThreadsX.map(x->AO(x), grid.points)
eval_AOs(grid::QuadGrid, AOs::Vector{AO}) = hcat(ThreadsX.map(Χμ->eval_AO(grid,Χμ), AOs)...)
