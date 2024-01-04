#!/usr/local/bin/julia

# John Eargle (mailto: jeargle at gmail.com)
# uveldtsim

using uveldt

using Dates

using ArgParse
using Printf


"""
    read_file_and_simulate(filename)

Simulate a system specified by a config file.

# Arguments
- filename: YAML config file
"""
function read_file_and_simulate(filename)

end


"""
    main()

Entrypoint for uveldt simulation script.
"""
function main()
    aps = ArgParseSettings()
    @add_arg_table! aps begin
        "configfile"
            help = "YAML system configuration file"
            required = true
    end

    parsed_args = parse_args(ARGS, aps)

    read_file_and_simulate(parsed_args["configfile"])
end

main()
