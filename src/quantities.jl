function quadripolar_moment(grid::QuadGrid, Ψs::Matrix{T}, R::T) where {T<:Real}
    ρ = sum(map(x->x .^2, eachcol(Ψs))) # density associated to Ψ

    quadmoment = zeros(3,3)
    quadmoment[end] = (R^2)/2 # nuclear contribution for the zz quadripolar moment

    for i in 1:3
        for j in 1:3
            δᵢⱼ = iszero(i-j) ? one(i) : zero(i)
            quadmoment_n_terms = map(enumerate(grid.points)) do (iX, X)
                grid.weights[iX]*( 3*X[i]*X[j] - δᵢⱼ*norm(X)^2 )*ρ[iX]
            end
            quadmoment[i, j] -= sum(quadmoment_n_terms)
        end
    end
    quadmoment
end
