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
eval_AO(grid::QuadGrid, AO::AO) = ThreadsX.map(x->AO(x), grid.points)
eval_AOs(grid::QuadGrid, AOs::Vector{AO}) = hcat(ThreadsX.map(Χμ->eval_AO(grid,Χμ), AOs)...)
