function run_helfem(Z1::String, Z2::String, Rbond::T;
                    multiplicity=1,
                    method="HF",
                    output_dir= Z1==Z2 ? Z1*"2_data" : Z1*Z2*"_data",
                    helfem_path,
                    write_hdf5=true,
                    helfem_kwargs...) where {T<:Real}
    # List of helfem routines
    helfem_commands = [joinpath(helfem_path, "objdir/src", command) for command in
                       ["diatomic_cbasis", "diatomic", "diatomic_dgrid"]]

    !isdir(output_dir) && mkdir(output_dir)

    # 1) Run cbasis routine to optimize the discretization basis to run helfem
    command = `$(helfem_commands[1]) --Z1=$(Z1) --Z2=$(Z2) --Rbond=$(Rbond) --angstrom=0`
    # run command and parse last line of the output
    cmd_output = split(read(command, String), "\n"; keepempty=false)[end]

    # 2) Run the .diatomic command
    # extract the helfem_kwargs as strings
    additional_kwargs = ["--$(key)=$(value)" for (key,value) in
                         zip(keys(helfem_kwargs), values(helfem_kwargs))]
    cmd_in = [helfem_commands[2], String.(split(cmd_output," "; keepempty=false))...,
              "--M=$(multiplicity)", "--method=$(method)", "--save=$(output_dir)/helfem_$(Rbond).chk",
              additional_kwargs...]
    cmd_output = read(Cmd(cmd_in), String)
    cmd_output = split(cmd_output,"\n",keepempty=false)

    # parse output to check kinetic and total energy
    parsed_output = parse_helfem_output(cmd_output)

    # 3) Extract grid data from checkfile
    if write_hdf5
        output_file = "$(output_dir)/helfem_$(Rbond).hdf5"
        cmd_in = [helfem_commands[3], "--load=$(output_dir)/helfem_$(Rbond).chk",
                  "--output=$(output_file)"]
        read(Cmd(cmd_in), String)

        # Clean the working dir
        rm("$(output_dir)/helfem_$(Rbond).chk")
        isfile("fort.9") && rm("fort.9")

        # Add parsed data to the HDF5 file.
        h5write(output_file, "Kinetic energy", parsed_output.e_kin)
        h5write(output_file, "Total energy", parsed_output.e_tot)
        h5write(output_file, "Electronic dipole", parsed_output.μ_elec)
        h5write(output_file, "Nuclear dipole", parsed_output.μ_nuc)
        h5write(output_file, "Total dipole", parsed_output.μ_tot)
        h5write(output_file, "Electronic quadrupole", parsed_output.Q_elec)
        h5write(output_file, "Nuclear quadrupole", parsed_output.Q_nuc)
        h5write(output_file, "Total quadrupole", parsed_output.Q_tot)
        h5write(output_file, "Hellmann-Feynman force", parsed_output.f_HF)
        return nothing
    else
        isfile("fort.9") && rm("fort.9")
        rm(output_dir; recursive=true)
        return parsed_output
    end
    nothing
end

function parse_helfem_output(cmd_output)
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

    # parse forces
    force_string = split(only( filter(x->contains(x, "Hellmann-Feynman"), cmd_output)), " ")[end]
    f_HF = parse(Float64, force_string)

    (; e_kin, e_tot, μ_elec, μ_nuc, μ_tot, Q_elec, Q_nuc, Q_tot, f_HF)
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
function reference_mo_coeffs(data::Dict{String, Any})
    ΨA = data["orba.re"] .+ im .* data["orba.im"]
    ΨB = data["orbb.re"] .+ im .* data["orbb.im"]
    (;α=ΨA, β=ΨB)
end

function reference_grad_mo_coeffs(data::Dict{String, Any})
    TΨA = data["Torba.re"] .+ im .* data["Torba.im"]
    TΨB = data["Torbb.re"] .+ im .* data["Torbb.im"]
    (;α=TΨA, β=TΨB)
end

function read_helfem_data(basis::String, datadir::String)
    # Extract raw data
    @assert(isdir(datadir))
    files = joinpath.(Ref(datadir), filter(x->!startswith(x, "_"), readdir(datadir)))
    for file in files
        @info "reading $(file)"
    end
    read_helfem_data(basis, files)
end
function read_helfem_data(basis::String, files::Vector{String})
    # Extract raw data
    output_data = (;)
    interatomic_distances = Float64[]
    reference_MOs = []
    reference_∇MOs = []
    grids = QuadGrid[]
    energies = Float64[]
    dipoles = Float64[]
    quadrupoles = Float64[]
    forces = Float64[]

    elements = nothing

    # Run through all JSON file.
    for filename in files # joinpath.(Ref(datadir), readdir(datadir))
        h5open(filename) do file
            data = read(file)
            # Extract Elements and grid a single time
            if (isempty(interatomic_distances))
                A, B = extract_elements(data, basis)
                elements = [A,B]
                output_data = merge(output_data, (;elements))
            end
            # Extracta and normalize reference eigenfunctions
            grid = QuadGrid(data)
            R = data["Rh"]*2

            # DEBUG: Only for closed-shell systems where α and β functions are the same (2 everywhere)
            Ψ = reference_mo_coeffs(data).α
            TΨ = reference_grad_mo_coeffs(data).α

            # Safely remove zero imaginary part
            @assert (norm(imag.(Ψ)) < 1e-10) && (norm(imag.(TΨ)) < 1e-10) "Non zero imaginary pat"
            Ψ = real.(Ψ)
            TΨ = real.(TΨ)

            # Check that the MOs are orthonormal and check kinetic term
            # precision. The tols are fixed for H2. Might break also check quadrupole.
            @assert norm(dot(grid, Ψ, Ψ) - I) < 1e-8 ""*
                "Reference eigenfunctions are not orthonormal"
            @assert norm(sum(diag(2*dot(grid, Ψ, TΨ))) - data["Kinetic energy"]) < 1e-6 ""*
                "Kinetic energies do not corresponds"
            test_quad = norm(quadrupole_moment(grid, elements..., Ψ, R, verbose=false)[end]
                             - data["Total quadrupole"])
            if test_quad > 1e-5
                @warn "Low quadrupole precision for interatomic distance $(R)!\n"*
                    "Distance to ref: $(test_quad)"
            end

            # Add data to reference dict
            push!(reference_MOs, Ψ)
            push!(reference_∇MOs, TΨ)
            push!(interatomic_distances, R)
            push!(grids, grid)
            push!(energies, data["Total energy"])
            ("Total dipole" ∈ keys(data)) && push!(dipoles, data["Total dipole"])
            push!(quadrupoles, data["Total quadrupole"])
            push!(forces, data["Hellmann-Feynman force"])
        end
    end

    # Return all data as a NamedTuple
    merge(output_data, (; interatomic_distances, reference_MOs, reference_∇MOs, grids,
                        basis, energies, dipoles, quadrupoles, forces))
end
