#!/usr/local/bin/julia --project=..

# John Eargle (mailto: jeargle at gmail.com)
# uveldtsim

using uveldt

using Dates

using ArgParse
using Printf


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

    setup_simulation(parsed_args["configfile"])
end

main()
