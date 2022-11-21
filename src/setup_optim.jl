"""
Define the inequality constraints in Optim.jl conventions.
"""
function setup_bounds(grid::QuadGrid{T}, A::Element{T}, B::Element{T};
                      num∫tol=1e-10, # maximum error of the numerical quadrature
                      ) where {T <: Real}
    # Maximal possible ζ such that ∫exp(-ζ*|x|²) is integrated with precsion "num∫tol".
    ζ_max = compute_spread_lim(grid; num∫tol)
    c_max = 1e3

    # Setup bounds
    nA, nB = length(vec(A)), length(vec(B))
    n_ζA, n_ζB = sum(A.shape_exps), sum(B.shape_exps)
    upper = vcat(ζ_max .* ones(n_ζA),      # upper bounds on spreads for A
                 c_max .* ones(nA-n_ζA),   # ------------ on ctr coeffs for A
                 ζ_max .* ones(n_ζB),      # ------------ on spreads for B
                 c_max .* ones(nB-n_ζB))   # ------------ on ctr coeffs for B     
    lower = vcat(zeros(n_ζA),              # lower bounds on spreads for the element A
                 -c_max .* ones(nA-n_ζA),  # ------------ on ctr coeffs for A
                 zeros(n_ζB),              # ------------ on spreads for B
                 -c_max .* ones(nB-n_ζB))  # ------------ on ctr coeffs for B
    upper, lower
end

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

    # Run through all JSON file. Only interatomic distance Rh changes for each file.
    for filename in joinpath.(Ref(datadir), readdir(datadir))
        h5open(filename) do file
            data = read(file)
            # Extract Elements and grid a single time
            if (isempty(Rhs))
                A, B = extract_elements(data, basis)
                Elements = [A,B]
                grid = QuadGrid(data)
                output_data = merge(output_data, (;Elements, grid))
            end
            
            push!(Rhs, data["Rh"])
            push!(Ψs_ref, reference_eigenvectors(data)[1])
        end       
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, basis))
end

function launch_optimization(ref_data; num∫tol=1e-7)
    # Inequality constraints of the basis optimization
    upper, lower = setup_bounds(ref_data.grid, ref_data.Elements...; num∫tol)
    # Print optimization parameters
    @info "GTO basis optimization\n"*
        "basis: $(ref_data.basis)\n"* "ζ_max: $(upper[1])\n"*
        "Maximum numerical integration error: $(num∫tol)"

    A, B = ref_data.Elements
    grid = ref_data.grid

    # Define objective function
    function j2opt(X::Vector{T}) where {T<:Real}
        accu = zero(T)
        for (i, Rh) in enumerate(ref_data.Rhs)
            accu += j_L2_diatomic(X, A, B, [0., 0., -Rh], [0., 0., Rh], ref_data.Ψs_ref[i], grid)
        end
        accu
    end

    # Start from standard basis
    X_init = vcat(vec(A), vec(B))
    if maximum(X_init) > upper[1]
        error("Some GTO cannot be integrated on the grid with"*
              " an error inferior to $(num∫tol). Lower num∫tol at your own risks.")
    end

    # Optimize with Fminbox routine (adapted for inequality constraints optimization)
    res = optimize(j2opt, lower, upper, X_init, Fminbox(LBFGS()),
                   Optim.Options(show_trace=true, extended_trace=true), autodiff=:forward);    
end
