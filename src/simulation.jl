# John Eargle (mailto: jeargle at gmail.com)
# uveldt.simulation

"""
Membrane-bound compartment containing Molecules and a Genome.  Chemical Reactions happen
within.
"""
struct Cell
    uuid::UUID
    genome::Genome
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    energy::Int64

    function Cell(genome::Genome)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(undef, 2)
        molecule_counts[1] = Dict{AbstractString, Int64}()
        molecule_counts[2] = Dict{AbstractString, Int64}()
        new(UUIDs.uuid4(), genome, molecule_counts, 0)
    end
end


"""
Unique location in a Veldt.  Contains two Molecule count buffers and possibly a single Cell.
VeldtPoint time evolution is controlled by a Veldt.
"""
mutable struct VeldtPoint
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    cell::Union{Cell, Nothing}

    function VeldtPoint(; molecules=[], cell=nothing)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(undef, 2)
        molecule_counts[1] = Dict{AbstractString, Int64}()
        molecule_counts[2] = Dict{AbstractString, Int64}()
        for mol in molecules
            molecule_counts[1][mol] = 0
            molecule_counts[2][mol] = 0
        end

        return new(molecule_counts, cell)
    end
end


"""
Multidimensional array representing locations that can hold Molecules and Cells.
"""
mutable struct Veldt
    dims::Array{Int64, 1}   # num points per dimension; 2 or 3 dimensions
    points  # Array with length(dims) dimensions holding VeldtPoints
    current # buffer with current (not future) data, 1 or 2
    molecule_counts::Dict{AbstractString, Int64}  # Molecules in VeldtPoints
    cell_molecule_counts::Dict{AbstractString, Int64}  # Molecules in Cells

    function Veldt(dims::Array{Int64, 1})
        if length(dims) == 2
            points = [[VeldtPoint() for j in 1:dims[2]]
                      for i in 1:dims[1]]
        elseif length(dims) == 3
            points = [[[VeldtPoint() for k in 1:dims[3]]
                       for j in 1:dims[2]]
                      for i in 1:dims[1]]
        else
        end
        new(dims, points, 1, Dict{AbstractString, Int64}(), Dict{AbstractString, Int64}())
    end
end


"""
    init_molecules(veldt, coord, molecule_counts)

Initialize the Molecule count for a specific location in a Veldt.

# Arguments
- `veldt::Veldt`:
- `coord::Array{Int64, 1}`:
- `molecule_counts::Dict{AbstractString, Int64}`:
"""
function init_molecules(veldt::Veldt, coord::Array{Int64, 1}, molecule_counts::Dict{String, Int64})
    for (mol, count) in molecule_counts
        if haskey(veldt.molecule_counts, mol)
            veldt.molecule_counts[mol] += count
        else
            veldt.molecule_counts[mol] = count
        end

        if length(coord) == 2
            vp = veldt.points[coord[1]][coord[2]]
            vp.molecule_counts[1][mol] = count
        elseif length(coord) == 3
            vp = veldt.points[coord[1]][coord[2]][coord[3]]
            vp.molecule_counts[1][mol] = count
        end
    end
end


"""
    add_cell(veldt, coord, cell)

Add a cell to a specific location in a Veldt.

# Arguments
- `veldt::Veldt`:
- `coord::Array{Int64, 1}`:
- `cell::Cell`:
"""
function add_cell(veldt::Veldt, coord::Array{Int64, 1}, cell::Cell)
    for (mol, count) in cell.molecule_counts[1]
        if haskey(veldt.cell_molecule_counts, mol)
            veldt.cell_molecule_counts[mol] += count
        else
            veldt.cell_molecule_counts[mol] = count
        end

        if length(coord) == 2
            vp = veldt.points[coord[1]][coord[2]]
            vp.cell = cell
        elseif length(coord) == 3
            vp = veldt.points[coord[1]][coord[2]][coord[3]]
            vp.cell = cell
        end
    end
end


"""
    get_neighbors(veldt, coord)

Return VeldtPoints neighboring VeldtPoint at coord.
2D: [x, -x, y, -y]
3D: [x, -x, y, -y, z, -z]
If neighbor does not exist, defaults to VeldtPoint at coord.

# Arguments
- `veldt::Veldt`:
- `coord::Array{Int64, 1}`:

# Returns
- `Array`: Array of VeldtPoints
"""
function get_neighbors(veldt::Veldt, coord::Array{Int64, 1})

    if length(veldt.dims) == 2
        x, y = coord
        neighbors = [veldt.points[x][y] for i in 1:4]

        if x > 1
            neighbors[1] = veldt.points[x-1][y]
        end

        if x < veldt.dims[1]
            neighbors[2] = veldt.points[x+1][y]
        end

        if y > 1
            neighbors[3] = veldt.points[x][y-1]
        end

        if y < veldt.dims[2]
            neighbors[4] = veldt.points[x][y+1]
        end
    elseif length(veldt.dims) == 3
        x, y, z = coord
        neighbors = [veldt.points[x][y][z] for i in 1:6]

        if x > 1
            neighbors[1] = veldt.points[x-1][y][z]
        end

        if x < veldt.dims[1]
            neighbors[2] = veldt.points[x+1][y][z]
        end

        if y > 1
            neighbors[3] = veldt.points[x][y-1][z]
        end

        if y < veldt.dims[2]
            neighbors[4] = veldt.points[x][y+1][z]
        end

        if z > 1
            neighbors[5] = veldt.points[x][y][z-1]
        end

        if z < veldt.dims[3]
            neighbors[6] = veldt.points[x][y][z+1]
        end
    end

    return neighbors
end


"""
    print_molecule_counts(veldt)

Print molecule counts per VeldtPoint in a Veldt.

# Arguments
- `veldt::Veldt`:
"""
function print_molecule_counts(veldt)
    current = veldt.current

    mol_count = 0
    for (mol, count) in veldt.molecule_counts
        mol_count += count
    end
    @printf "total: %d\n" mol_count

    if length(veldt.dims) == 2
        for i in 1:veldt.dims[1]
            @printf "  %d:" i
            for j in 1:veldt.dims[2]
                vp = veldt.points[i][j]
                mol_count = 0
                for (mol, count) in vp.molecule_counts[current]
                    # @printf "  %s: %d\n" mol count
                    mol_count += count
                end
                # @printf "  (%d,%d): %d" i j mol_count
                @printf "%4d" mol_count
            end
            println()
        end
    elseif length(veldt.dims) == 3
        for i in 1:veldt.dims[1]
            @printf "  slice %d\n" i
            for j in 1:veldt.dims[2]
                @printf "  %d:" j
                for k in 1:veldt.dims[3]
                    vp = veldt.points[i][j][k]
                    mol_count = 0
                    for (mol, count) in vp.molecule_counts[current]
                        mol_count += count
                    end
                    # @printf "  (%d,%d,%d): %d\n" i j k mol_count
                    @printf "%4d" mol_count
                end
                println()
            end
            println()
        end
    end

end


"""
    setup_veldt(filepath)

Create a Veldt from a YAML setup file.

# Arguments
- `filepath`: path to YAML setup file

# Returns
- `Veldt`:
"""
function setup_veldt(filepath)
    # Move to setup file directory.
    dir = dirname(filepath)
    base_dir = pwd()
    if dir != ""
        cd(dir)
    end

    filename = basename(filepath)
    setup = YAML.load(open(filename))

    # Build Veldt.
    if haskey(setup, "dimensions")
        dimensions = setup["dimensions"]
    end

    veldt = Veldt(dimensions)

    # Build Chemistry.
    if haskey(setup, "chemistry")
        chemistry_file = setup["chemistry"]
    end
    chemistry = setup_chemistry(chemistry_file)


    # Build Genomes.
    genomes = Dict{String, Genome}()
    if haskey(setup, "genomes")
        genome_files = setup["genomes"]
        for genome_file in genome_files
            genome_info = read_fasta(genome_file)
            for (name, uuid, genome_str) in genome_info
                if uuid == nothing
                    @printf "  genome: %s\n" name
                    genomes[name] = Genome(name, genome_str, chemistry)
                else
                    @printf "  genome: %s %s\n" name uuid
                    genomes[name] = Genome(uuid, name, genome_str, chemistry)
                end
            end
        end
    end

    @printf "  genomes: %s\n" genomes

    # Build Molecules.
    if haskey(setup, "molecules")
        for mol_info in setup["molecules"]
            if haskey(mol_info, "name")
                mol_name = mol_info["name"]
            end

            if haskey(mol_info, "locations")
                for loc_info in mol_info["locations"]
                    mol_loc = loc_info["location"]
                    mol_count = loc_info["count"]
                    init_molecules(veldt, mol_loc, Dict(mol_name => mol_count))
                end
            end
        end
    end

    # Build Cells.
    if haskey(setup, "cells")
        for cell_info in setup["cells"]
            if haskey(cell_info, "genome")
                genome_name = cell_info["genome"]
            end
            cell = Cell(genomes[genome_name])

            if haskey(cell_info, "molecules")
                for mol_info in cell_info["molecules"]
                    mol_name = mol_info["name"]
                    mol_count = mol_info["count"]
                    cell.molecule_counts[1][mol_name] = mol_count
                end
            end

            if haskey(cell_info, "location")
                location = cell_info["location"]
                add_cell(veldt, location, cell)
            end
        end
    end

    if dir != ""
        cd(base_dir)
    end

    return veldt
end


"""
    step2D(veldt::Veldt)

Simulate one timestep for a 2D Veldt.

# Arguments
- `veldt::Veldt`:

# Returns
- `Veldt`: Veldt after one timestep
"""
function step2D(veldt::Veldt)

    current = veldt.current
    next = (current == 1) ? 2 : 1

    # println("  current: ", current)

    # VeldtPoints
    # for each VeldtPoint
    for i in 1:veldt.dims[1]
        for j in 1:veldt.dims[2]
            vp = veldt.points[i][j]
            # Choose molecules to move in each direction (transport or diffusion)
            #   2D: 4 directions; +1 into Cell
            neighbors = get_neighbors(veldt, [i, j])
            for (mol, count) in vp.molecule_counts[current]
                # Move molecules to neighboring VeldtPoints or internal Cell; next buffer
                for k in 1:count
                    chosen = rand([1, 2, 3, 4])
                    if haskey(neighbors[chosen].molecule_counts[next], mol)
                        neighbors[chosen].molecule_counts[next][mol] += 1
                    else
                        neighbors[chosen].molecule_counts[next][mol] = 1
                    end
                end
                vp.molecule_counts[current][mol] = 0
            end

            if vp.cell != nothing

                # Choose molecules to react
                # Run reactions
                # Put products in next buffer
                # Choose molecules to move (transport or diffusion)
                #   1 direction out of Cell
                # Move molecules to containing VeldtPoint; next buffer
            end
        end
    end

    # Switch current buffer
    veldt.current = next

    return veldt
end


"""
    step3D(veldt::Veldt)

Simulate one timestep for a 3D Veldt.

# Arguments
- `veldt::Veldt`:

# Returns
- `Veldt`: Veldt after one timestep
"""
function step3D(veldt::Veldt)

    current = veldt.current
    next = (current == 1) ? 2 : 1

    # println("  current: ", current)

    # VeldtPoints
    # for each VeldtPoint
    for i in 1:veldt.dims[1]
        for j in 1:veldt.dims[2]
            for k in 1:veldt.dims[3]
                vp = veldt.points[i][j][k]
                # Choose molecules to move in each direction (transport or diffusion)
                #   3D: 6 directions; +1 into Cell
                neighbors = get_neighbors(veldt, [i, j, k])
                for (mol, count) in vp.molecule_counts[current]
                    # Move molecules to neighboring VeldtPoints or internal Cell; next buffer
                    for k in 1:count
                        chosen = rand([1, 2, 3, 4, 5, 6])
                        if haskey(neighbors[chosen].molecule_counts[next], mol)
                            neighbors[chosen].molecule_counts[next][mol] += 1
                        else
                            neighbors[chosen].molecule_counts[next][mol] = 1
                        end
                    end
                    vp.molecule_counts[current][mol] = 0
                end

                #   if Cell
                if vp.cell != nothing
                    #     Choose molecules to react
                    #     Run reactions
                    #     Put products in next buffer
                    #     Choose molecules to move (transport or diffusion)
                    #       1 direction out of Cell
                    #     Move molecules to containing VeldtPoint; next buffer
                end
            end
        end
    end

    # Switch current buffer
    veldt.current = next

    return veldt
end


"""
    simulate(veldt::Veldt, numsteps)

Simulate one or more timesteps for a Veldt.

# Arguments
- `veldt::Veldt`: Veldt to simulate
- `numsteps`: number of timesteps to simulate

# Returns
- `Veldt`: Veldt after one or more timesteps
"""
function simulate(veldt::Veldt, numsteps)

    if length(veldt.dims) == 2
        print_molecule_counts(veldt)
        println()
        for i in 1:numsteps
            println("step: ", i)
            step2D(veldt)
            print_molecule_counts(veldt)
            println()
        end
    elseif length(veldt.dims) == 3
        print_molecule_counts(veldt)
        println()
        for i in 1:numsteps
            println("step: ", i)
            step3D(veldt)
            print_molecule_counts(veldt)
            println()
        end
    end

    # Record state.

    return veldt
end


"""
    setup_simulation(filepath)

Load and run a simulation from a YAML setup file.

# Arguments
- `filepath`: path to YAML setup file
"""
function setup_simulation(filepath)
    # Move to setup file directory.
    dir = dirname(filepath)
    base_dir = pwd()
    if dir != ""
        cd(dir)
    end

    filename = basename(filepath)
    setup = YAML.load(open(filename))

    # Build Veldt.
    if haskey(setup, "veldt")
        veldt_file = setup["veldt"]
    end
    veldt = setup_veldt(veldt_file)

    # Simulation parameters
    step_count = 0
    if haskey(setup, "steps")
        step_count = setup["steps"]
    end

    if dir != ""
        cd(base_dir)
    end

    simulate(veldt, step_count)
end
