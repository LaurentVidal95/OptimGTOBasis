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

    # parse output to check kinetic and total energy
    parsed_energies = filter(x->(contains(x, "Total") || contains(x, "Kinetic")) &&
                             (contains(x, "energy")), cmd_output)[end-1:end]
    e_kin = parse(Float64, split(parsed_energies[1], " ")[end])
    e_tot = parse(Float64, split(parsed_energies[2], " ")[end])

    # 3) Extract grid data from checkfile
    output_file = "$(output_dir)/helfem_$(Rbond).hdf5"
    cmd_in = [helfem_commands[3], "--load=$(output_dir)/helfem_$(Rbond).chk",
              "--output=$(output_file)"]
    read(Cmd(cmd_in), String)

    # Clean the working dir
    # rm("$(output_dir)/helfem_$(Rbond).chk")
    isfile("fort.9") && rm("fort.9")

    # Add e_tot to the HDF5 file.
    h5write(output_file, "Kinetic energy", e_kin)
    h5write(output_file, "Total energy", e_tot)
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
