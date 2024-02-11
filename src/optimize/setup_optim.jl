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
    files = joinpath.(Ref(datadir), filter(x->!startswith(x, "_"), readdir(datadir)))
    for file in files
        @info "reading $(file)"
    end
    extract_ref_data(basis, files)
end
function extract_ref_data(basis::String, files::Vector{String})
    # Extract raw data
    output_data = (;)
    Rhs = []
    Ψs_ref = []
    grids = QuadGrid[]
    Energies = Float64[]
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
            # Extracta and normalize reference eigenfunctions
            grid = QuadGrid(data)
            Ψs = reference_eigenvectors(data)[1]
            @assert(norm(imag.(Ψs)) < 1e-10) # Check that the functions are real
            grid_norm = sqrt.(diag(dot(grid, Ψs, Ψs)))
            Ψs = real.((1 ./ grid_norm)' .* Ψs)

            push!(Rhs, data["Rh"])
            push!(grids, grid)
            push!(Ψs_ref, Ψs)
            push!(Energies, data["Total energy"])
        end
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, grids, basis, Energies))
end

function default_starting_point(ref_data)
    A, B = ref_data.Elements
    X_start = A==B ? vec(A) : vcat(vec(A), vec(B)) # Initial guess
    X_start
end

"""
Define the inequality constraints in Optim.jl conventions.
"""
function setup_bounds!(model::Model, ref_data, ζ_max::T;
                       exps_tol=1e-2,
                       X_start=default_starting_point(ref_data)) where {T <: Real}
    c_max = 500  # Defaut maximum absolute value of ctr coefficient

    # Extract needed parameters from ref_data
    A, B = ref_data.Elements
    nA, nB = length(vec(A)), length(vec(B))
    n_params = A==B ? nA : nA+nB
    @assert length(X_start)==n_params
    n_ζA, n_ζB = sum(A.shape_exps), sum(B.shape_exps)


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

function maximum_spread(ref_data, criterion::OptimizationCriterion; default_max_spread=200)
    isa(criterion, EnergyCriterion) && (return default_max_spread)
    findmin([compute_spread_lim(grid, Rh; criterion.gridtol) for (grid, Rh) in
             zip(criterion.grids, (1/2) .* criterion.interatomic_distances)])[1]
end

function setup_optim_model(ref_data, criterion::OptimizationCriterion;
                           X_start=default_starting_point(ref_data))
    @assert Bool(sum(isa.(Ref(criterion), [ProjectionCriterion, EnergyCriterion]))) ""*
        "Choose between energy and density criterion"

    criterion_type = typeof(criterion)==ProjectionCriterion ? :projection : :energy
    @info "Setup model for the $(criterion_type) criterion"
    model = Model(Ipopt.Optimizer)

    # Setup constrained variables
    ζ_max = maximum_spread(ref_data, criterion)
    @info "Maximum spread: $(ζ_max)"
    setup_bounds!(model, ref_data, ζ_max; X_start)

    # Extract elements and handle the A=B case.
    A, B = ref_data.Elements
    n_params = A==B ? length(vec(A)) : length(vec(A)) + length(vec(B))

    j2opt(X::T...) where {T<:Real} = objective_function(criterion, A, B, X...)
    if isa(criterion, ProjectionCriterion)
        register(model, :j2opt, n_params, j2opt; autodiff=true)
    else
        ∇j2opt!(∇J, X::T...) where {T<:Real} = grad_objective_function!(criterion, A, B, ∇J, X...)
        register(model, :j2opt, n_params, j2opt, ∇j2opt!)
    end
    X = model[:X]
    @NLobjective(model, Min, j2opt(X...))
    model
end

# TODO
function extract_optimized_data(model::Model)
    nothing
end
