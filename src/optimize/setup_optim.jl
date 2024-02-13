# Beware some imaginary parts are non negligeable.
function reference_eigenvectors(data::Dict{String, Any})
    ΨA = data["orba.re"] .+ im .* data["orba.im"]
    ΨB = data["orbb.re"] .+ im .* data["orbb.im"]
    ΨA, ΨB
end

function reference_kinetic(data::Dict{String, Any})
    TΨA = data["Torba.re"] .+ im .* data["Torba.im"]
    TΨB = data["Torbb.re"] .+ im .* data["Torbb.im"]
    TΨA, TΨB
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
    Rhs = Float64[]
    Ψs_ref = []
    TΨs_ref = []
    grids = QuadGrid[]
    Energies = Float64[]

    # Run through all JSON file.
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
            # DEBUG: Only alpha here....
            Ψs = reference_eigenvectors(data)[1]
            TΨs = reference_kinetic(data)[1]
            # Sanity checks
            # Check that the functions are real
            @assert(norm(imag.(Ψs)) < 1e-10) 
            @assert(norm(imag.(TΨs)) < 1e-10)
            Ψs = real.(Ψs)
            TΨs = real.(TΨs)
            # Check that Ψs are orthonormal and check kinetic term precision
            @assert norm(dot(grid, Ψs, Ψs) - I) < 1e-8 ""*
                "Reference eigenfunctions are not orthonormal"
            @assert norm(sum(2*dot(grid, Ψs, TΨs)) - data["Kinetic energy"]) < 1e-6

            # Add data to reference dict
            push!(TΨs_ref, TΨs)
            push!(Rhs, data["Rh"])
            push!(grids, grid)
            push!(Ψs_ref, Ψs)
            push!(Energies, data["Total energy"])
        end
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, TΨs_ref, grids, basis, Energies))
end

function set_starting_point(ref_data; guess=:bse)
    @assert guess ∈ (:bse, :random)
    A, B = ref_data.Elements
    n_params = A==B ? length(vec(A)) : sum(length.(vec(A), vec(B)))
    if guess==:bse
        (A==B) && return vec(A)
        return vcat(vec(A), vec(B))
    else
        return rand(n_params)
    end
    error("Bug in starting point")
end


function launch_Optim(ref_data, criterion::OptimizationCriterion, X_guess;
                      maxiter=500, solver=ConjugateGradient())
    A, B = ref_data.Elements
    f(X) = objective_function(criterion, A, B, X...)
    g!(∇E, X) = grad_objective_function!(criterion, A, B, ∇E, X...)
    # Choose between energy and projection criterion
    if false #isa(criterion, EnergyCriterion)
        return optimize(f, g!, X_guess, solver,
                        Optim.Options(show_trace=true, iterations=maxiter))
    else
        return optimize(f, X_guess, solver, Optim.Options(show_trace=true, iterations=maxiter))
    end
    error("Not supposed to happen")
end

function launch_Ipopt(ref_data, criterion::OptimizationCriterion, X_guess;
                      maxiter=50)
    model = ipopt_setup_optim_model(ref_data, criterion, X_guess)
    set_optimizer_attribute(model, "line_search_method", "filter")
    set_optimizer_attribute(model, "max_iter", maxiter)
    optimize!(model)
    model
end

# TODO
function extract_optimized_data(model::Model)
    nothing
end
