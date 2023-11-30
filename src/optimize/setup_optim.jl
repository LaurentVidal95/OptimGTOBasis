using JuMP
using Ipopt

# Beware some imaginary parts are non negligeable.
function reference_eigenvectors(data::Dict{String, Any})
    ΨA = data["orba.re"] .+ im .* data["orba.im"]
    ΨB = data["orbb.re"] .+ im .* data["orbb.im"]
    ΨA, ΨB
end

function extract_ref_data(basis::String, datadir::String)
    # Extract raw data
    @assert(isdir(datadir))
    files = joinpath.(Ref(datadir), readdir(datadir))
    extract_ref_data(basis, files)
end
function extract_ref_data(basis::String, files::Vector{String})
    # Extract raw data
    output_data = (;)
    Rhs = []
    Ψs_ref = []
    grids = QuadGrid[]
    Energies = []
    normalize_col(tab) = hcat(normalize.(eachcol(tab))...)

    # Run through all JSON file. Only interatomic distance Rh changes for each file.
    for filename in files # joinpath.(Ref(datadir), readdir(datadir))
        h5open(filename) do file
            data = read(file)
            # Extract Elements and grid a single time
            if (isempty(Rhs))
                A, B = extract_elements(data, basis)
                Elements = [A,B]
                output_data = merge(output_data, (;Elements))
            end            
            push!(Rhs, data["Rh"])
            push!(Ψs_ref, normalize_col(reference_eigenvectors(data)[1]))
            push!(grids, QuadGrid(data))
            push!(Energies, data["Total energy"])
        end       
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, grids, basis, Energies))
end

"""
Define the inequality constraints in Optim.jl conventions.
"""
function setup_bounds!(model::Model, ref_data, ζ_max::T;
                       exps_tol=1e-2) where {T <: Real}
    c_max = 500  # Defaut maximum absolute value of ctr coefficient

    # Extract needed parameters from ref_data
    A, B = ref_data.Elements
    X_start = A==B ? vec(A) : vcat(vec(A), vec(B)) # Initial guess
    nA, nB = length(vec(A)), length(vec(B))
    n_ζA, n_ζB = sum(A.shape_exps), sum(B.shape_exps)

    n_params = A==B ? nA : nA+nB

    # Setup variables of the problem and initial guess.
    @variable(model, X[i=1:n_params], start=X_start[i])
    # Defined boxed constraints with JuMP conventions
    @constraint(model, [i=1:n_ζA],  exps_tol ≤     X[i]    ≤ ζ_max, base_name="spread_A")
    @constraint(model, [i=1:nA-n_ζA], -c_max ≤  X[i+n_ζA]  ≤ c_max, base_name="ctr_A")
    if !(A==B)
        @constraint(model, [i=1:n_ζB],  exps_tol ≤   X[nA+i]   ≤ ζ_max, base_name="spread_B")
        @constraint(model, [i=1:nB-n_ζB], -c_max ≤ X[nA+n_ζB+i] ≤ c_max, base_name="ctr_B")
    end
    nothing
end

"""
Setup the optimization problem using JuMP framework.
"""
function setup_optim_model(ref_data; num∫tol=1e-7)
    model = Model(Ipopt.Optimizer)

    # Compute maximum spread and setup constrained variables
    ζ_max = findmin([compute_spread_lim(grid, Rh; num∫tol)
                     for (grid, Rh) in zip(ref_data.grids, ref_data.Rhs)])[1]
    @info "Maximum possible spread for given integration grids: $(ζ_max)"
    setup_bounds!(model, ref_data, ζ_max)

    # Extract elements and handle the A=B case.
    A, B = ref_data.Elements
    n_params = A==B ? length(vec(A)) : length(vec(A)) + length(vec(B))
    
    # Define and register objective function 
    # function j2opt(X::T...) where {T<:Real}
    #     Y = collect(X)
    #     accu = zero(T)
    #     for (i, Rh) in enumerate(ref_data.Rhs)
    #         # accu += j_L2_diatomic(Y, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Ψs_ref[i],
    #         #                       ref_data.grids[i])
    #         accu += j_E_diatomic(Y, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Energies[i])
    #     end
    #     accu
    # end
    function j2opt(X::T...) where {T<:Real}
        Y = collect(X)
        sum([j_E_diatomic(Y, A, B, [0., 0., -Rh], [0., 0., Rh], E)
                         for (Rh, E) in zip(ref_data.Rhs, ref_data.Energies)])
    end
    function ∇j2opt!(∇J, X::T...) where {T<:Real}
        Y = collect(X)
        ∇Y = zero(Y)
        for (i, Rh) in enumerate(ref_data.Rhs)
            ∇Y .+= ∇j_E_diatomic(Y, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Energies[i])
        end
        for i in 1:length(Y)
            ∇J[i] = ∇Y[i]
        end
        return nothing
    end
    register(model, :j2opt, n_params, j2opt, ∇j2opt!) #; autodiff=true)
    X = model[:X]
    @NLobjective(model, Min, j2opt(X...))
    model
end

# TODO
function extract_optimized_data(model::Model)
    nothing
end
