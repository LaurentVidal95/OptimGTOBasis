function run_helfem(Z1::String, Z2::String, Rbond::T;
                    multiplicity=1,
                    method="HF",
                    output_dir= Z1==Z2 ? Z1*"2_data" : Z1*Z2*"_data",
                    verbose=true,
                    helfem_path) where {T<:Real}
    # List of helfem routines
    helfem_commands = [joinpath(helfem_path, "objdir/src", command) for command in
                       ["diatomic_cbasis", "diatomic", "diatomic_dgrid"]]

    !isdir(output_dir) && mkdir(output_dir)
    run_fct = verbose ? cmd->run(cmd) : cmd->read(cmd, String)
    
    # 1) Run cbasis routine to optimize the discretization basis to run helfem
    command = `$(helfem_commands[1]) --Z1=$(Z1) --Z2=$(Z2) --Rbond=$(Rbond) --angstrom=0`
    # run command and parse last line of the output
    cmd_output = split(read(command, String), "\n"; keepempty=false)[end]
    
    # 2) Run the .diatomic command
    cmd_in = [helfem_commands[2], String.(split(cmd_output," "; keepempty=false))...,
              "--M=$(multiplicity)", "--method=$(method)", "--save=$(output_dir)/helfem_$(Rbond).chk"]
    run_fct(Cmd(cmd_in))

    # 3) Extract data from checkfile
    cmd_in = [helfem_commands[3], "--load=$(output_dir)/helfem_$(Rbond).chk",
              "--output=$(output_dir)/helfem_$(Rbond).hdf5"]
    run_fct(Cmd(cmd_in))

    # Clean the working dir
    rm("$(output_dir)/helfem_$(Rbond).chk")
    isfile("fort.9") && rm("fort.9")
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
