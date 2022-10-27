"""
Contains the exponents [ζ_i] and contracting coeffs [c_i] such that
`` R(X) = ∑_i c_i exp(-ζi*|X|^2)``
"""
struct RadialPart{T<:Real}
    exps   ::Vector{T}
    coeffs ::Vector{T}
end
(Rnl::RadialPart)(X) = Rnl.coeffs'exp.(map(ζ -> exp(-ζ*norm(X)^2), Rnl.exps))

"""
Ψ_{nlm} = Y_lm * R_nl
• The radial part R_nl is defined above
• The real spherical harmonics is defined with the SphericalHarmonicExpansions
  package.
"""
struct AO{T<:Real, F}
    Rnl    ::RadialPart{T}
    Ylm    ::F
    center ::Vector{T}
end
function AO(l::Int64, m::Int64, exps::Vector{T}, coeffs::Vector{T},
            center::Vector{T}) where {T<:Real}
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
function construct_AOs(X::Element;
                       position=zeros(T, 3), # [0., 0., ± Rh]
                       mmax) where {T<:Real}
    AO_basis = AO[]
    @info "Added AOs: "
    # Run other n, l and m ∈ {-l,...,l} and add each corresponding AO
    for (i_shell, shell) in enumerate(X.shells)
        l = i_shell-1
        for n in 1:size(shell.coeffs,2)
            for m in -min(mmax, l):min(mmax,l)
                label = ["s","p","d","f","g"][i_shell]
                println(@sprintf("%-10s %-6s",  "Χ_{$(n+l)$l$m}", "$(n+l)$label"))
                push!(AO_basis, AO(l, m, shell.exps, vec(shell.coeffs[:,n]), position))
            end
        end
    end
    AO_basis
end
eval_AOs(AOs::Vector{AO}, grid::QuadGrid) = hcat(map(Χμ->Χμ.(grid.points), AOs)...)


