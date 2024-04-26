function eval_abs_ms(mol::PyObject)
    # Compute AO labels for molecule X₂
    all_labels = mol.ao_labels(fmt=false)
    N_aos = div(length(all_labels),2)
    labels = all_labels[1:N_aos] # AO label on first X

    # Compute abs of the associated quantum number m
    multipopfirst!(tab, N) = [popfirst!(tab) for x in tab[1:N]]
    abs_ms = Int[]
    while !isempty(labels)
        l = split(labels[1][3],"")[end]
        n_aos_nl=Dict("s"=>1, "p"=>3, "d"=>5, "f"=>7, "g"=>9,
                      "h"=>11)
        aos_nl = multipopfirst!(labels, n_aos_nl[l])
        if l=="s"
            append!(abs_ms, [0])
        elseif l=="p"
            append!(abs_ms, [1,1,0])
        else # l=d,f,g,...
            # Compute numbers of aos with same n and l
            half_n_ao_nl = div(length(aos_nl), 2)
            append!(abs_ms, abs.(collect(-half_n_ao_nl:half_n_ao_nl)))
        end
    end
    vcat(abs_ms, abs_ms)
end

function eval_AOs(grid::QuadGrid, A::Element{T1}, B::Element{T1},
                  R::T2; deriv=0) where {T1, T2<:Real}
    @assert deriv ∈ (0,2) # never use first derivative
    basis = Dict([A.name => basis_string([A]),
                  B.name => basis_string([B])])
    tmp_mol = pyscf.M(;atom="$(A.name) 0.0 0.0 -$(R/2);
                             $(B.name) 0.0 0.0 $(R/2)",
                      basis,
                      unit="bohr"
                      )
    X = pyscf.dft.numint.eval_ao(tmp_mol, grid.points; deriv)
    if deriv==2
        # X is a tensor. The first coordinates are
        # 1: AOs on the grid
        # 2-4: Derivative along x, y and z axes
        # 5-10: Second derivatives along xx, xy, xz, yy, yz, zz, hence the 5, 8, and last.
        # The laplacian is sum(xx, yy, zz)
        X = X[5, :, :] .+ X[8,:,:] .+ X[10,:,:]
    end

    # filter X with mmax: X[:,i]=0 if |m(X[:,i])| > mmax
    id_selected_aos = [i for (i,m) in enumerate(eval_abs_ms(tmp_mol)) if m≤grid.mmax]
    X[:,id_selected_aos]
end

function overlap(grid::QuadGrid, A::Element{T}, B::Element{T}, R; norm_type=:L²) where {T<:Real}
    basis = Dict([X.name => basis_string([X]) for X in (A,B)])
    tmp_mol = pyscf.M(;atom="$(A.name) 0.0 0.0 -$(R/2);
                             $(B.name) 0.0 0.0 $(R/2)",
                      basis,
                      unit="bohr"
                      )
    id_selected_aos = [i for (i,m) in enumerate(eval_abs_ms(tmp_mol)) if m≤grid.mmax]
    S = tmp_mol.intor("int1e_ovlp")
    if norm_type==:L²
        return S[id_selected_aos, id_selected_aos]
    end
    kin = 2*tmp_mol.intor("int1e_kin")
    (S + kin)[id_selected_aos, id_selected_aos]
end
