using JuMP, Ipopt

# Beware some imaginary parts are non negligeable.
function reference_eigenvectors(data::Dict{String, Any})
    ΨA = data["orba.re"] .+ im .* data["orba.im"]
    ΨB = data["orbb.re"] .+ im .* data["orbb.im"]
    ΨA, ΨB
end

function extract_ref_data(basis::String, datadir::String)
    # Extract raw data
    @assert(isdir(datadir))
    output_data = (;)
    Rhs = []
    Ψs_ref = []
    grids = QuadGrid[]
    # Run through all JSON file. Only interatomic distance Rh changes for each file.
    for filename in joinpath.(Ref(datadir), readdir(datadir))
        h5open(filename) do file
            data = read(file)
            # Extract Elements and grid a single time
            if (isempty(Rhs))
                A, B = extract_elements(data, basis)
                Elements = [A,B]
                output_data = merge(output_data, (;Elements))
            end            
            push!(Rhs, data["Rh"])
            push!(Ψs_ref, reference_eigenvectors(data)[1])
            push!(grids, QuadGrid(data))
        end       
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, grids, basis))
end

"""
Define the inequality constraints in Optim.jl conventions.
"""
function setup_bounds!(model::Model, ref_data, ζ_max::T) where {T <: Real}
    c_max = 500  # Defaut maximum ctr coefficient

    # Extract needed parameters from ref_data
    A, B = ref_data.Elements
    X_start = vcat(vec(A), vec(B))
    nA, nB = length(vec(A)), length(vec(B))
    n_ζA, n_ζB = sum(A.shape_exps), sum(B.shape_exps)

    # Setup variables of the problem with constraints in JuMP conventions
    @variable(model, X[i=1:nA+nB], start=X_start[i])
    @constraint(model, [i=1:n_ζA],    0≤ X[i]      ≤ ζ_max, base_name="spread_A")
    @constraint(model, [i=1:nA-n_ζA], 0≤ X[i+n_ζA] ≤ c_max, base_name="ctr_A")
    @constraint(model, [i=1:n_ζB],    0≤ X[nA+i]      ≤ ζ_max, base_name="spread_B")
    @constraint(model, [i=1:nB-n_ζB], 0≤ X[nA+n_ζB+i] ≤ c_max, base_name="ctr_B")
    nothing
end

function setup_optim_model(ref_data; num∫tol=1e-7)
    model = Model(Ipopt.Optimizer)

    # Compute maximum spread and setup constrained variables
    ζ_max = findmin([compute_spread_lim(grid, Rh; num∫tol)
                     for (grid, Rh) in zip(ref_data.grids, ref_data.Rhs)])[1]
    setup_bounds!(model, ref_data, ζ_max)

    # Extract elements
    A, B = ref_data.Elements
    n_params = length(vec(A)) + length(vec(B))

    # Define and register objective function 
    function j2opt(X::T...) where {T<:Real}
        # @show [x.value for x in X]
        Y = collect(X)
        accu = zero(T)      
        for (i, Rh) in enumerate(ref_data.Rhs)
            accu += j_L2_diatomic(Y, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Ψs_ref[i],
                                  ref_data.grids[i])
        end
        accu
    end
    register(model, :j2opt, n_params, j2opt; autodiff=true)
    @NLobjective(model, Min, j2opt(model[:X]...))
    model
end

# function launch_optimization(ref_data; num∫tol=1e-7,
#                              method=LBFGS,
#                              kwargs...)
#     # Inequality constraints of the basis optimization
#     i_rough_grid = findmin([compute_spread_lim(grid; num∫tol) for grid in ref_data.grids])[2]
#     upper, lower = setup_bounds(ref_data.grids[i_rough_grid], ref_data.Rhs[i],
#                                 ref_data.Elements...; num∫tol)
#     # Print optimization parameters
#     # @info "GTO basis optimization\n"*
#     #     "basis: $(ref_data.basis)\n"* "ζ_max: $(upper[1])\n"*
#     #     "Maximum numerical integration error: $(num∫tol)"

#     A, B = ref_data.Elements

#     # Define objective function
#     function j2opt(X::Vector{T}) where {T<:Real}
#         # @show [x.value for x in X]
#         accu = zero(T)
#         for (i, Rh) in enumerate(ref_data.Rhs)
#             accu += j_L2_diatomic(X, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Ψs_ref[i],
#                                   ref_data.grids[i])
#         end
#         accu
#     end

#     # Start from standard basis
#     X_init = vcat(vec(A), vec(B))
#     if maximum(X_init) > upper[1]
#         error("Some GTO cannot be integrated on the grid with"*
#               " an error inferior to $(num∫tol). Lower num∫tol at your own risks.")
#     end
#     @show lower, upper
#     # Optimize with Fminbox routine (adapted for inequality constraints optimization)
#     res = optimize(j2opt, lower, upper, X_init, Fminbox(method(;kwargs...)),
#                    Optim.Options(show_trace=true, extended_trace=true, iterations=20),
#                    autodiff=:forward);    
# end
