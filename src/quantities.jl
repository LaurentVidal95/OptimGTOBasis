function quadrupole(grid::QuadGrid, A::Element, B::Element, Ψs::Matrix{T}, R::T) where {T<:Real}
    # DEBUG: only for closed shell... Otherwise you have to track the alpha and beta part..
    ρ = 2*sum(map(x->x .^2, eachcol(Ψs))) # density associated to Ψ

    quadmoment = zeros(3,3)
    # DEBUG: to be implemented (Z1 + Z2)*(R^2)/4
    quadmoment[end] = (A.charge + B.charge)*(R^2)/4 # nuclear contribution for the zz quadripolar moment

    for i in 1:3
        for j in 1:3
            δᵢⱼ = iszero(i-j) ? one(i) : zero(i)
            quadmoment_n_terms = map(enumerate(grid.points)) do (iX, X)
                (1/2)*grid.weights[iX]*( 3*X[i]*X[j] - δᵢⱼ*norm(X)^2 )*ρ[iX]
            end
            quadmoment[i, j] -= sum(quadmoment_n_terms)
        end
    end
    quadmoment
end
