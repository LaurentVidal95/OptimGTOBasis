function quadrupole(grid::QuadGrid, A::Element, B::Element, Ψs::Matrix{T}, R::T) where {T<:Real}
    # DEBUG: only for closed shell... Otherwise you have to track the alpha and beta part..
    ρ = 2*sum(map(x->x .^2, eachcol(Ψs))) # density associated to Ψ

    quadmoment = zeros(3,3)
    # nuclear contribution for the zz quadripolar moment
    quadmoment[end] = (A.charge + B.charge)*(R^2)/4

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

function compare_quadrupoles_homoatomic(file::String, bases::Vector{String}, bases_names, ref_data)
    @warn "Only for closed shell homoatomic molecule"
    @assert isfile(file)
    data = open(JSON3.read, file)

    output = Dict{Any, Any}()
    # run on bases
    for (basis, name) in zip(bases, bases_names)
        A = Element(Vector(data[name]), ref_data.Elements[1])
        P = ProjectionCriterion(ref_data)

        output[name] = zeros(length(bases), 2)
        for i in 1:length(ref_data.Ψs_ref)

            # Extract ref data for the i-th file
            R = P.interatomic_distances[i]
            grid = P.grids[i]
            Ψ = P.reference_functions[i]
            TΨ = P.reference_kinetics[i]

            # Compute orthogonal projections in L² and H¹ norms
            Ψ_proj_L2 = orthogonal_projection(A, A, R, Ψ, TΨ, grid; norm_type=:L²)
            # Ψ_proj_H1 = orthogonal_projection(A, A, R, Ψ, TΨ, grid; norm_type=:H¹)

            # Compute corresponding approximate quadrupoles Qzz
            quad_L2 = quadrupole(grid, A, A, Ψ_proj_L2, R)[end]
            quad_H1 = quadrupole(grid, A, A, Ψ_proj_H1, R)[end]
            output[name][i,1] = quad_L2
            output[name][i,2] = quad_H1
        end
    end
    output
end
