function set_starting_point(ref_data; guess=:bse)
    @assert guess ∈ (:bse, :random)
    A, B = ref_data.elements
    n_params = A==B ? length(vec(A)) : sum(length.(vec(A), vec(B)))
    if guess==:bse
        if (A==B)
            X_guess = vec(A)
            n_exps = len_exps(A)
            X_guess[1:n_exps] = log.(X_guess[1:n_exps])
            return X_guess
        else
            n_exps = len_exps.([A,B])
            XA_guess = vec(A)
            XB_guess = vec(B)
            XA_guess[1:n_exps[1]] .= log.(XA_guess[1:n_exps[1]])
            XB_guess[1:n_exps[1]] .= log.(XB_guess[1:n_exps[1]])
            vcat(XA_guess, XB_guess)
        end
    else
        return rand(n_params)
    end
    error("Bug in starting point")
end


function launch_Optim(ref_data, criterion::OptimizationCriterion, X_guess;
                      maxiter=500, solver=LBFGS(), tol=5e-9)
    A, B = ref_data.elements
    f(X) = objective_function(criterion, A, B, X...; kwargs...)
    # Choose between energy and projection criterion
    return optimize(f, X_guess, solver, Optim.Options(show_trace=true, iterations=maxiter,
                                                      g_tol=tol))
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
