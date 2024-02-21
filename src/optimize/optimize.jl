function optimize_AO_basis(basis::String, datadir::String, criterion_type;
                           optimizer=:Optim, # Choose between Optim.jl and JuMP (Ipopt)
                           maxiter=100,
                           norm_type=:L²,
                           guess=:bse,
                           kwargs...)
    ref_data = read_helfem_data(basis, datadir)
    criterion = criterion_type(ref_data; norm_type)
    X_guess = set_starting_point(ref_data; guess)

    if (optimizer==:Ipopt)
        return launch_Ipopt(ref_data, criterion, X_guess; maxiter, kwargs...)
    else
        return launch_Optim(ref_data, criterion, X_guess; maxiter, kwargs...)
    end
    error("Not supposed to happen")
end
