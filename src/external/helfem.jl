function run_helfem(Z1::String, Z2::String, Rbond::T;
                    multiplicity=1,
                    method="HF",
                    output_dir= Z1==Z2 ? Z1*"2_data" : Z1*Z2*"_data",
                    helfem_path) where {T<:Real}
    # List of helfem routines
    helfem_commands = [joinpath(helfem_path, "objdir/src", command) for command in
                       ["diatomic_cbasis", "diatomic", "diatomic_dgrid"]]

    !isdir(output_dir) && mkdir(output_dir)

    # 1) Run cbasis routine to optimize the discretization basis to run helfem
    command = `$(helfem_commands[1]) --Z1=$(Z1) --Z2=$(Z2) --Rbond=$(Rbond) --angstrom=0`
    # run command and parse last line of the output
    cmd_output = split(read(command, String), "\n"; keepempty=false)[end]

    # 2) Run the .diatomic command
    cmd_in = [helfem_commands[2], String.(split(cmd_output," "; keepempty=false))...,
              "--M=$(multiplicity)", "--method=$(method)", "--save=$(output_dir)/helfem_$(Rbond).chk"]
    cmd_output = read(Cmd(cmd_in), String)
    cmd_output = split(cmd_output,"\n",keepempty=false)
    println(cmd_output)

    # parse output to check kinetic and total energy
    parsed_energies = filter(x->(contains(x, "Total") || contains(x, "Kinetic")) &&
                             (contains(x, "energy")), cmd_output)[end-1:end]
    e_kin = parse(Float64, split(parsed_energies[1], " ")[end])
    e_tot = parse(Float64, split(parsed_energies[2], " ")[end])

    # parse dipole moments
    parsed_dipole = split.(filter(x->contains(x, "dipole"), cmd_output), " "; keepempty=false)
    μ_elec = parse(Float64, parsed_dipole[1][4])
    μ_nuc = parse(Float64, parsed_dipole[2][4])
    μ_tot = parse(Float64, parsed_dipole[3][4])

    # pase quadrupole moments
    parsed_quadrupole = split.(filter(x->contains(x,"quadrupole"), cmd_output),
                               " "; keepempty=false)
    Q_elec = parse(Float64, parsed_quadrupole[1][4])
    Q_nuc = parse(Float64, parsed_quadrupole[2][4])
    Q_tot = parse(Float64, parsed_quadrupole[3][4])


    # 3) Extract grid data from checkfile
    output_file = "$(output_dir)/helfem_$(Rbond).hdf5"
    cmd_in = [helfem_commands[3], "--load=$(output_dir)/helfem_$(Rbond).chk",
              "--output=$(output_file)"]
    read(Cmd(cmd_in), String)

    # Clean the working dir
    rm("$(output_dir)/helfem_$(Rbond).chk")
    isfile("fort.9") && rm("fort.9")

    # Add energies and quadrupole to the HDF5 file.
    h5write(output_file, "Kinetic energy", e_kin)
    h5write(output_file, "Total energy", e_tot)
    h5write(output_file, "Electronic dipole", μ_elec)
    h5write(output_file, "Nuclear dipole", μ_nuc)
    h5write(output_file, "Total dipole", μ_tot)
    h5write(output_file, "Electronic quadrupole", Q_elec)
    h5write(output_file, "Nuclear quadrupole", Q_nuc)
    h5write(output_file, "Total quadrupole", Q_tot)

    nothing
end

function generate_reference_data(Z1::String, Z2::String, bond_lengths::Vector{T};
                                 kwargs...) where {T<:Real}
    NR = length(bond_lengths)
    for (i,R) in enumerate(bond_lengths)
        @info "Reference computation: $i/$NR"
        run_helfem(Z1, Z2, R; kwargs...)
    end
    nothing
end

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
    quadrupoles = Float64[]
    Elements = nothing

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
            Rh = data["Rh"]

            # DEBUG: Only for closed-shell systems (2 everywhere)
            Ψs = reference_eigenvectors(data)[1]
            TΨs = reference_kinetic(data)[1]
            @assert(norm(imag.(Ψs)) < 1e-10)
            @assert(norm(imag.(TΨs)) < 1e-10)
            Ψs = real.(Ψs)
            TΨs = real.(TΨs)

            # Check that Ψs are orthonormal and check kinetic term precision
            # The tols are fixed for H2. Might break
            # also check quadrupole
            @assert norm(dot(grid, Ψs, Ψs) - I) < 1e-8 ""*
                "Reference eigenfunctions are not orthonormal"
            @assert norm(sum(diag(2*dot(grid, Ψs, TΨs))) - data["Kinetic energy"]) < 1e-6 ""*
                "Kinetic energies do not corresponds"
            test_quad = norm(quadrupole(grid, Elements..., Ψs, Rh*2)[end] - data["Total quadrupole"])
            if test_quad > 1e-5
                @warn "Low quadrupole precision for interatomic distance $(Rh*2)!\n"*
                    "Distance to ref: $(test_quad)"
            end

            # Add data to reference dict
            push!(TΨs_ref, TΨs)
            push!(Rhs, Rh)
            push!(grids, grid)
            push!(Ψs_ref, Ψs)
            push!(Energies, data["Total energy"])
            push!(quadrupoles, data["Total quadrupole"])
        end
    end

    # Return all data as a NamedTuple
    merge(output_data, (; Rhs, Ψs_ref, TΨs_ref, grids, basis, Energies, quadrupoles))
end
