# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt

module uveldt

# Chemistry
include("chemistry.jl")
export Element, ElementTable, Bond, BondTable, Chemistry, Molecule
export ReactionType, Reaction

export get_bond, mass, setup_chemistry


export Gene, Genome, Cell, VeldtPoint, Veldt, init_molecules, add_cell
export parse_reaction, transcribe_gene, genome_string, find_genes
export add_snps, add_insertions, remove_deletions, cross_over
export is_pseudogene
export read_fasta, write_fasta
export setup_veldt

using Distributions
using Printf
using Random
using YAML


"""
Genes are strings bracketed by start '(' and stop ')' characters that code for specific
chemical reactions, pores, or transporters.

For chemical reactions, the internals consist of Molecule strings with join '*' and split
'/' operators specifying how the set of Molecules will be processed.  All operators act at
the same time.  Each join creates a new bond, and each split breaks an existing bond.  To
create a reaction string for the inverse of a given reaction, just switch all joins to
splits and vice versa.

# Examples
  (A*B): A + B -> AB
  (A/B): AB -> A + B
  (AB*C): AB + C -> ABC
  (AB/C): ABC -> AB + C
  (A*B/C): A + BC -> AB + C
  (A*B*C): A + B + C -> ABC
  (A*A*A): A + A + A -> AAA

Pores are specified by single Molecule strings with no join or split operators.  If a pore
is active, Molecules of that type are allowed to diffuse freely into and out of the cell.

# Examples
  (A): A pore
  (AB): AB pore
  (AAA): AAA pore

Transporters are similar to pores except that they use energy to move Molecules across the
cell membrane.  A transporter string must start with a join or split operator followed by
a single Molecule string.

# Examples
  (*A): A(out) -> A(in)
  (/A): A(in) -> A(out)
  (*AB): AB(out) -> AB(in)
  (/AAA): AAA(in) -> AAA(out)

"""
struct Gene
    location::Int64
    string::AbstractString
    transcript::AbstractString
    chemistry::Chemistry

    function Gene(location::Int64, gene_string::AbstractString, chemistry::Chemistry)
        transcript = transcribe_gene(gene_string, chemistry)
        new(location, gene_string, transcript, chemistry)
    end
end


"""
String specifying all genes in the cell.  Genes are regions of the Genome that fit the Gene
pattern.  A Genome can have a lot of intergenic space that doesn't code for anything.
"""
struct Genome
    name::AbstractString
    string::AbstractString
    chemistry::Chemistry
    genes::Array{Gene, 1}

    Genome(name::AbstractString, string::AbstractString, chemistry::Chemistry) = new(name, string, chemistry)

    function Genome(name::AbstractString, size::Int64, chemistry::Chemistry)
        Genome(name, genome_string(size, chemistry), chemistry)
    end
end


"""
Membrane-bound compartment containing Molecules and a Genome.  Chemical Reactions happen
within.
"""
struct Cell
    genome::Genome
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    energy::Int64

    function Cell(genome::Genome)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(undef, 2)
        molecule_counts[1] = Dict{AbstractString, Int64}()
        molecule_counts[2] = Dict{AbstractString, Int64}()
        new(genome, molecule_counts, 0)
    end
end


"""
Unique location in a Veldt.  Contains two Molecule count buffers and possibly a single Cell.
VeldtPoint time evolution is controlled by a Veldt.
"""
mutable struct VeldtPoint
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    cell::Cell

    # function VeldtPoint()
    #     molecule_counts = Array{Dict{AbstractString, Int64}, 1}(2)
    #     molecule_counts[1] = Dict{AbstractString, Int64}()
    #     molecule_counts[2] = Dict{AbstractString, Int64}()
    #     new(molecule_counts)
    # end

    function VeldtPoint(; molecules=[], cell=nothing)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(undef, 2)
        molecule_counts[1] = Dict{AbstractString, Int64}()
        molecule_counts[2] = Dict{AbstractString, Int64}()
        for mol in molecules
            molecule_counts[1][mol] = 0
            molecule_counts[2][mol] = 0
        end
        if cell == nothing
            return new(molecule_counts)
        else
            return new(molecule_counts, cell)
        end
    end
end


"""
Multidimensional array representing locations that can hold Molecules and Cells.
"""
struct Veldt
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
    init_molecules(veldt, coord, molecules_counts)

Initialize the Molecule count for a specific location in a Veldt.

# Arguments
- veldt::Veldt
- coord::Array{Int64, 1}
- molecule_counts::Dict{AbstractString, Int64}
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
- veldt::Veldt
- coord::Array{Int64, 1}
- cell::Cell
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
    parse_reaction(reaction_string)

Parse the reactants and products from a reaction string.

# Arguments
- reaction_string::AbstractString

# Returns
- (reaction_type::ReactionType, reactants::Array, products::Array)
"""
function parse_reaction(reaction_string::AbstractString)
    if isnothing(match(r"[\*/]", reaction_string))
        reaction_type = pore
        reactants = [reaction_string]
        products = []
    elseif reaction_string[1] == '*'
        reaction_type = transport_in
        reactants = [replace(reaction_string, r"\*|/" => "")]
        products = []
    elseif reaction_string[1] == '/'
        reaction_type = transport_out
        reactants = [replace(reaction_string, r"\*|/" => "")]
        products = []
    else
        reaction_type = reaction
        reactants = split(replace(reaction_string, "/" => ""), "*")
        products = split(replace(reaction_string, "*" => ""), "/")
    end

    return (reaction_type, reactants, products)
end


"""
    transcribe_gene(gene_string, chemistry)

Convert a Gene string into a valid Reaction string by removing redundant directive
characters.

# Arguments
- gene_string::AbstractString
- chemistry::Chemistry

# Returns
- Reaction string
"""
function transcribe_gene(gene_string::AbstractString, chemistry::Chemistry)
    # modes for parsing
    #   0: in Molecule
    #   1: in directive
    mode = 0
    gene_chars = []
    elements = join(keys(chemistry.element_table.elements))

    if !startswith(gene_string, "(")
        println("Error: Gene must start with \"(\"")
    end

    if !endswith(gene_string, ")")
        println("Error: Gene must end with \")\"")
    end

    for el in gene_string[2:end-1]
        if mode == 0
            if el == '*' || el == '/'
                push!(gene_chars, el)
                mode = 1
            elseif !in(el, elements)
                push!(gene_chars, el)
                mode = 0
            end
        elseif mode == 1
            if !in(el, elements)
                push!(gene_chars, el)
                mode = 0
            elseif el != '*' && el != '/'
                @printf "Error: bad character in Gene string - %s" el
            end
        end
    end

    # Do not end string with a directive character.
    if mode == 1
        pop!(gene_chars)
    end

    return join(gene_chars)
end


"""
    genome_string(size, chemistry)

Randomly generate a valid Genome string.

# Arguments
- size::Int64
- chemistry::Chemistry

# Returns
- Genome string
"""
function genome_string(size::Int64, chemistry::Chemistry)
    genome = randstring(alphabet_string(chemistry), size)
    # return join(genome)
    return genome
end


"""
    find_genes(genome)

Search through a Genome for patterns that match possible Gene strings and build Genes from
them.

# Arguments
- genome::Genome

# Returns
- Array of Genes
"""
function find_genes(genome::Genome)
    rx = r"\([^\(]*?\)"
    gene_match = eachmatch(rx, genome.string)
    genes = [Gene(gm.offset, gm.match, genome.chemistry)
             for gm in gene_match]
    return genes
end


"""
    add_snps(genome, rate)

Add SNPs to Genome.

# Arguments
- genome::Genome
- rate

# Returns
- Mutated Genome String
"""
function add_snps(genome::Genome, rate)
    geom_dist = Geometric(rate)
    location = 1 + rand(geom_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        snp = string(rand(alphabet))
        @printf "  %d SNP %s->%s\n" location genome.string[location] snp
        push!(fragments, snp)
        start = location + 1
        location += rand(geom_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    add_insertions(genome, rate; size_param)

Add small insertions to Genome.

# Arguments
- genome::Genome
- rate
- size_param

# Returns
- Mutated Genome String
"""
function add_insertions(genome::Genome, rate; size_param=0.5)
    location_dist = Geometric(rate)
    location = 1 + rand(location_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        insert = randstring(alphabet, size)
        @printf "  %d insert %s\n" location insert
        push!(fragments, insert)
        start = location
        location += rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    remove_deletions(genome, rate)

Remove small deletions from Genome.

# Arguments
- genome::Genome
- rate
- size_param

# Returns
- Mutated Genome String
"""
function remove_deletions(genome::Genome, rate; size_param=0.5)
    geom_dist = Geometric(rate)
    location = 1 + rand(geom_dist)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        @printf "  %d delete %s\n" location genome.string[location:location+size]
        start = location + size
        location += rand(geom_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    add_duplications(genome, rate; size_param)

Add large duplications to Genome.

# Arguments
- genome::Genome
- rate
- size_param

# Returns
- Mutated Genome String
"""
function add_duplications(genome::Genome, rate; size_param=0.5)
    location_dist = Geometric(rate)
    location = 1 + rand(location_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        # insert = randstring(alphabet, size)
        # @printf "  %d insert %s\n" location insert
        # push!(fragments, insert)
        start = location
        location += rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    add_inversions(genome, rate; size_param)

Add large inversions to Genome.

# Arguments
- genome::Genome
- rate
- size_param

# Returns
- Mutated Genome String
"""
function add_inversions(genome::Genome, rate; size_param=0.5)
    location_dist = Geometric(rate)
    location = 1 + rand(location_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        # insert = randstring(alphabet, size)
        # @printf "  %d insert %s\n" location insert
        # push!(fragments, insert)
        start = location
        location += rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    add_translocations(genome, rate; size_param)

Add large translocations to Genome.

# Arguments
- genome::Genome
- rate
- size_param

# Returns
- Mutated Genome String
"""
function add_translocations(genome::Genome, rate; size_param=0.5)
    location_dist = Geometric(rate)
    location = 1 + rand(location_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        # insert = randstring(alphabet, size)
        # @printf "  %d insert %s\n" location insert
        # push!(fragments, insert)
        start = location
        location += rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    return join(fragments)
end


"""
    cross_over(genome1, genome2)

Cross over two Genomes to produce two child Genomes.

# Arguments
- genome1::Genome
- genome2::Genome

# Returns
- Tuple of mutated Genome Strings, (String, String)
"""
function cross_over(genome1::Genome, genome2::Genome)
    if genome1.chemistry != genome2.chemistry
        @printf "Error: chemistries must match between genome1 and genome2"
    end

    location1 = rand(0:length(genome1.string))
    location2 = rand(0:length(genome2.string))
    @printf "  location1: %d\n" location1
    @printf "  location2: %d\n" location2

    head_str = genome1.string[1:location1]
    tail_str = genome2.string[location2+1:end]
    @printf "  head: %s\n" head_str
    @printf "  tail: %s\n" tail_str
    child_str1 = head_str * tail_str

    head_str = genome2.string[1:location2]
    tail_str = genome1.string[location1+1:end]
    @printf "  head: %s\n" head_str
    @printf "  tail: %s\n" tail_str
    child_str2 = head_str * tail_str

    return (child_str1, child_str2)
end


"""
    is_pseudogene(gene)

Identify pseudogene.

# Arguments
- gene::Gene

# Returns
- Bool
"""
function is_pseudogene(gene::Gene)
    elements = join(keys(gene.chemistry.element_table.elements))

    if length(gene.transcript) == 0
        return true
    end

    m = match(r"([\*/]*[$(elements)])", gene.transcript)
    if m === nothing
        return true
    end

    return false
end


"""
    read_fasta(filename)

Read a FASTA file and return the information as tuples of names and a genome strings.

# Arguments
- filename

# Returns
- Array of tuples (name string, genome string)
"""
function read_fasta(filename)
    genome_info = []
    name = ""
    genome_str = ""
    first_genome = true
    open(filename, "r") do f
        for line in eachline(f)
            println(line)
            if line[1] == '>'
                if !first_genome
                    push!(genome_info, (name, genome_str))
                else
                    first_genome = false
                end
                name = line[2:end]
                genome_str = ""
            else
                genome_str *= line
            end
        end
    end
    push!(genome_info, (name, genome_str))
    return genome_info
end


"""
    write_fasta(genomes, filename)

Write a multi-FASTA file with the full genome strings for an array of Genomes.

# Arguments
- genomes
- filename
"""
function write_fasta(genomes, filename)
    f = open(filename, "w")
    line_length = 80
    for genome in genomes
        println(f, ">", genome.name)
        full_line_count = div(length(genome.string), line_length)  # integer division
        for i in 1:full_line_count
            println(f, genome.string[(i-1)*line_length+1:i*line_length])
        end
        println(f, genome.string[full_line_count*line_length+1:length(genome.string)])
    end
    close(f)
end


"""
    write_fasta(genome, filename)

Write a FASTA file with the full genome string.

# Arguments
- genome::Genomes
- filename
"""
function write_fasta(genome::Genome, filename)
    write_fasta([genome], filename)
end




"""
    setup_veldt(filename)

Create a Veldt from a YAML setup file.

# Arguments
- filename: name of YAML setup file

# Returns
- Veldt
"""
function setup_veldt(filename)
    setup = YAML.load(open(filename))

    # build Veldt
    if haskey(setup, "dimensions")
        dimensions = setup["dimensions"]
    end

    veldt = Veldt(dimensions)

    # build Chemistry
    if haskey(setup, "chemistry")
        chemistry_file = setup["chemistry"]
    end
    chemistry = setup_chemistry(chemistry_file)


    # build Genomes
    genomes = Dict{String, Genome}()
    if haskey(setup, "genomes")
        genome_files = setup["genomes"]
        for genome_file in genome_files
            genome_info = read_fasta(genome_file)
            for (name, genome_str) in genome_info
                @printf "  genome: %s\n" name
                genomes[name] = Genome(name, genome_str, chemistry)
            end
        end
    end

    @printf "  genomes: %s\n" genomes

    # build Cells
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

    return Veldt
end


end
