# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt.evolution

"""
Node in a Phylogeny.
"""
struct PhyloNode
    genomeName::AbstractString
    genomeUuid::UUID

    function PhyloNode(genome)
        new(genome)
    end
end


"""
Edge in a Phylogeny.
"""
struct PhyloEdge
    parent::PhyloNode
    child::PhyloNode
    mutations

    function PhyloNode(genome)
        new(genome)
    end
end


"""
Phylogenic graph for multiple, related Genomes.
"""
struct Phylogeny
    nodes::Array{PhyloNode, 1}
    edges::Array{PhyloEdge, 1}
    roots::Array{PhyloNode, 1}
    leaves::Array{PhyloNode, 1}

    function Phylogeny(nodes, edges)
        new(nodes, edges)
    end
end


"""
Selection parameters for evolution algorithm.

Convert string to function name:
  julia> s = Symbol("sin")
  :sin

  julia> f = getfield(Main, s)
  sin (generic function with 10 methods)

  julia> f(3)
  0.1411200080598672

  julia> sin(3)
  0.1411200080598672
"""
struct SelectionParams
    fitness_function
    select_count::Int64
    select_fraction::Float64
    fitness_threshold::Float64

    function SelectionParams(fitness_function; select_count=0, select_fraction=0.0,
                             fitness_threshold=0.0)
        new(fitness_function, select_count, select_fraction, fitness_threshold)
    end
end


"""
Substitution matrix for SNVs.
"""
struct SubstitutionMatrix
    alphabet::AbstractString
    substitutions::Dict{Char, Categorical}

    function SubstitutionMatrix(alphabet, substitutions)
        new(alphabet, substitutions)
    end
end


"""
Mutation parameters for evolution algorithm.
"""
struct MutationParams
    snv_rate::Float64
    substitution_matrix::Union{SubstitutionMatrix, Nothing}
    insertion_rate::Float64
    insertion_size::Float64
    deletion_rate::Float64
    deletion_size::Float64
    duplication_rate::Float64
    duplication_size::Float64
    inversion_rate::Float64
    inversion_size::Float64
    translocation_rate::Float64
    translocation_size::Float64
    crossing_over::Bool

    function MutationParams(; snv_rate=0.0, substitution_matrix=nothing,
                            insertion_rate=0.0, insertion_size=0.5,
                            deletion_rate=0.0, deletion_size=0.5,
                            duplication_rate=0.0, duplication_size=0.02,
                            inversion_rate=0.0, inversion_size=0.02,
                            translocation_rate=0.0, translocation_size=0.02,
                            crossing_over=false)
        new(snv_rate, substitution_matrix, insertion_rate,
            insertion_size, deletion_rate, deletion_size,
            duplication_rate, duplication_size, inversion_rate,
            inversion_size, translocation_rate, translocation_size, crossing_over)
    end
end


"""
    select_genomes(genomes, params)

Take an Array of Genomes and return a subset of them based on
calculated fitness scores.

# Arguments
- genomes::Array{Genome, 1}
- params::SelectionParams

# Returns
- Array{Genome, 1} sub array of selected Genomes
"""
function select_genomes(genomes, params::SelectionParams)
    uuid_to_genome = Dict(genome.uuid => genome
                          for genome in genomes)

    # Calculate fitness scores.
    fitnesses = [(params.fitness_function(genome), genome.uuid)
                 for genome in genomes]

    # Sort by fitness.
    sort!(fitnesses, by = x -> x[1], rev=true)

    # Select final Genomes.
    if params.select_count > 0
        # by count
        if params.select_count < length(fitnesses)
            fitnesses = fitnesses[1:params.select_count]
        end
    elseif params.select_fraction > 0.0 && params.select_fraction <= 1.0
        # by fraction
        select_count = Int(params.select_fraction * length(fitnesses))
        if select_count < length(fitnesses)
            fitnesses = fitnesses[1:select_count]
        end
    elseif params.fitness_threshold > 0.0
        # by fitness threshold
        filter!(fitness -> fitness >= params.fitness_threshold, fitnesses)
    end

    selected_genomes = [uuid_to_genome[uuid]
                        for (fitness, uuid) in fitnesses]

    return selected_genomes
end


"""
    select_cells(cells, params)

Take an Array of Cells and return a subset of them based on
calculated fitness scores.

# Arguments
- cells::Array{Cell, 1}
- params::SelectionParams

# Returns
- Array{Cell, 1} sub array of selected Cells
"""
function select_cells(cells, params::SelectionParams)
    uuid_to_cell = Dict(cell.uuid => cell
                        for cell in cells)

    # Calculate fitness scores.
    fitnesses = [(params.fitness_function(cell), cell.uuid)
                 for cell in cells]

    # Sort by fitness.
    sort!(fitnesses, by = x -> x[1], rev=true)

    # Select final Cells.
    if params.select_count > 0
        # by count
        if params.select_count < length(fitnesses)
            fitnesses = fitnesses[1:params.select_count]
        end
    elseif params.select_fraction > 0.0 && params.select_fraction <= 1.0
        # by fraction
        select_count = Int(params.select_fraction * length(fitnesses))
        if select_count < length(fitnesses)
            fitnesses = fitnesses[1:select_count]
        end
    elseif params.fitness_threshold > 0.0
        # by fitness threshold
        filter!(fitness -> fitness >= params.fitness_threshold, fitnesses)
    end

    selected_cells = [uuid_to_cell[uuid]
                      for (fitness, cell) in fitnesses]

    return selected_cells
end


#####
# Fitness functions
#####

function gene_count(genome)
    return length(find_genes(genome))
end

function genome_length(genome)
    return length(genome.string)
end



"""
    read_substition_matrix(filename)

Read a substitution matrix file and return a SubstitutionMatrix.

Substitution matrix files have the form:

     A    B    C    (    )    *    /
A    0    5    5    1    1    1    1
B    5    0    5    1    1    1    1
C    5    5    0    1    1    1    1
(    2    2    2    0    1    1    1
)    2    2    2    1    0    1    1
*    2    2    2    1    1    0    1
/    2    2    2    1    1    1    0

where 'A', 'B', and 'C' are Elements in the Chemistry, and '(', ')',
'*', and '/' are the standard operator characters.  The numbers Mij
give odds that a character in row i will turn into a character in
column j.  SNV target selection already takes into account which bases
will remain the same so the matrix diagonal (Mii) should normally be
set to 0.

# Arguments
- filename

# Returns
- SubstitutionMatrix
"""
function read_substitution_matrix(filename, chemistry)
    alphabet = alphabet_string(chemistry)
    substitutions = Dict{Char, Categorical}()

    lines = readlines(filename)
    ordered_chars = join(split(lines[1]))
    for oc in ordered_chars
        if !(oc in alphabet)
            error("substitution characters must be in the provided Chemistry")
        end
    end

    for line in lines[2:end]
        line_contents = split(line)
        current_char = line_contents[1][1]
        # Use L1 norm to generate probabilities.
        odds = normalize([parse(Int64, lc) for lc in line_contents[2:end]], 1)
        substitutions[current_char] = Categorical(odds)
    end

    return SubstitutionMatrix(alphabet, substitutions)
end


"""
    substitute(sub_mat, original_char)

Pick a substitution character based on a SubstitutionMatrix.

# Arguments
- sub_mat::SubstitutionMatrix
- original_char::Char

# Returns
- Substitution Char
"""
function substitution(sub_mat::SubstitutionMatrix, original_char::Char)
    return sub_mat.alphabet[random(sub_mat.substitutions[original_char])]
end


"""
    add_snvs(genome, rate; sub_mat)

Add SNVs to Genome.

# Arguments
- genome::Genome
- rate
- sub_mat::SubstitutionMatrix

# Returns
- Mutated Genome
"""
function add_snvs(genome::Genome, rate; sub_mat::Union{SubstitutionMatrix, Nothing}=nothing)
    geom_dist = Geometric(rate)
    location = 1 + rand(geom_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        if sub_mat == nothing
            snv = string(rand(alphabet))
        else
            snv = string(substitute(sub_mat, genome.string[location]))
        end
        # @printf "  %d SNV %s->%s\n" location genome.string[location] snv
        push!(fragments, snv)
        start = location + 1
        location += rand(geom_dist)
    end

    push!(fragments, genome.string[start:end])

    return Genome(join(fragments), genome.chemistry)
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
        # @printf "  %d insert %s\n" location insert
        push!(fragments, insert)
        start = location
        location += rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    return Genome(join(fragments), genome.chemistry)
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
        # @printf "  %d delete %s\n" location genome.string[location:location+size]
        start = location + size
        location += rand(geom_dist)
    end

    push!(fragments, genome.string[start:end])

    return Genome(join(fragments), genome.chemistry)
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

    size_dist = Geometric(size_param)

    # Find and duplicate segments.
    duplications = []
    while location <= length(genome.string)
        size = 1 + rand(size_dist)
        if length(genome.string) < location+size
            size = length(genome.string) - location
        end
        duplication = genome.string[location:location+size]
        push!(duplications, duplication)
        # @printf "  %d duplicate %s\n" location duplication
        location += size + rand(location_dist)
    end

    # Find insertion points.
    insertion_points = []
    location = 1 + rand(location_dist)
    for duplication in duplications
        while location > length(genome.string)
            location -= length(genome.string)
        end
        push!(insertion_points, location)
        location += rand(location_dist)
    end

    # Make sure insertion_points are in order, but keep them
    # associated with their corresponding duplications.
    ordered_inserts = collect(zip(insertion_points, duplications))
    sort!(ordered_inserts, by = x -> x[1])

    # Insert duplications.
    start = 1
    fragments = []
    # @printf "  length(genome.string) %d\n" length(genome.string)
    for (insertion_point, duplication) in ordered_inserts
        # @printf "  start: %d insertion_point-1 %d duplication %s\n" start insertion_point-1 duplication
        push!(fragments, genome.string[start:insertion_point-1])
        push!(fragments, duplication)
        start = insertion_point
    end

    push!(fragments, genome.string[start:end])

    return Genome(join(fragments), genome.chemistry)
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
    geom_dist = Geometric(rate)
    location = 1 + rand(geom_dist)
    fragments = []
    start = 1

    size_dist = Geometric(size_param)

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        size = 1 + rand(size_dist)
        if length(genome.string) < location+size
            size = length(genome.string) - location
        end
        # @printf "  %d invert %s\n" location genome.string[location:location+size]
        push!(fragments, reverse(genome.string[location:location+size]))
        start = location + size
        location += rand(geom_dist)
    end

    push!(fragments, genome.string[start:end])

    return Genome(join(fragments), genome.chemistry)
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

    # @printf "  length(genome.string) %d\n" length(genome.string)

    # Find, duplicate, and remove segments.
    duplications = []
    while location <= length(genome.string)
        size = 1 + rand(size_dist)
        if length(genome.string) < location+size
            size = length(genome.string) - location
        end
        push!(fragments, genome.string[start:location-1])
        duplication = genome.string[location:location+size]
        push!(duplications, duplication)
        # @printf "  %d duplicate %s %d\n" location duplication length(duplication)
        start = location + size + 1
        location += size + rand(location_dist)
    end

    push!(fragments, genome.string[start:end])

    # Build intermediate genome string with segments removed.
    genome_string = join(fragments)

    # Find insertion points.
    insertion_points = []
    location = 1 + rand(location_dist)
    for duplication in duplications
        while location > length(genome_string)
            location -= length(genome_string)
        end
        push!(insertion_points, location)
        location += rand(location_dist)
    end

    # Make sure insertion_points are in order, but keep them
    # associated with their corresponding duplications.
    ordered_inserts = collect(zip(insertion_points, duplications))
    sort!(ordered_inserts, by = x -> x[1])

    # Insert duplications.
    start = 1
    fragments = []
    # @printf "  length(genome_string) %d\n" length(genome_string)
    for (insertion_point, duplication) in ordered_inserts
        # @printf "  start: %d insertion_point-1 %d duplication %s\n" start insertion_point-1 duplication
        push!(fragments, genome_string[start:insertion_point-1])
        push!(fragments, duplication)
        start = insertion_point
    end

    push!(fragments, genome_string[start:end])

    # @printf "  length(genome.string) %d\n" length(join(fragments))

    return Genome(join(fragments), genome.chemistry)
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
    # @printf "  location1: %d\n" location1
    # @printf "  location2: %d\n" location2

    head_str = genome1.string[1:location1]
    tail_str = genome2.string[location2+1:end]
    # @printf "  head: %s\n" head_str
    # @printf "  tail: %s\n" tail_str
    child_genome1 = Genome(head_str * tail_str, genome1.chemistry)

    head_str = genome2.string[1:location2]
    tail_str = genome1.string[location1+1:end]
    # @printf "  head: %s\n" head_str
    # @printf "  tail: %s\n" tail_str
    child_genome2 = Genome(head_str * tail_str, genome1.chemistry)

    return (child_genome1, child_genome2)
end


"""
    mutate(genomes, params)

Take an Array of cloned parent Genomes and create a set of child Genomes.

Replication is done prior to this step so the number of input Genomes
should equal the number of mutant child Genomes.

# Arguments
- genomes::Array{Genome, 1}

# Returns
- Array{Genome, 1} containing new child Genomes
"""
function mutate(genomes, params::MutationParams)
    child_genomes = []

    for genome in genomes
        if params.snv_rate > 0
            genome = add_snvs(genome, params.snv_rate;
                              sub_mat=params.substitution_matrix)
        end

        if params.insertion_rate > 0
            genome = add_insertions(genome, params.insertion_rate;
                                    size_param=params.insertion_size)
        end

        if params.deletion_rate > 0
            genome = remove_deletions(genome, params.deletion_rate;
                                      size_param=params.deletion_size)
        end

        if params.duplication_rate > 0
            genome = add_duplications(genome, params.duplication_rate;
                                      size_param=params.duplication_size)
        end

        if params.inversion_rate > 0
            genome = add_inversions(genome, params.inversion_rate;
                                    size_param=params.inversion_size)
        end

        if params.translocation_rate > 0
            genome = add_translocations(genome, params.translocation_rate;
                                        size_param=params.translocation_size)
        end

        push!(child_genomes, genome)
    end

    # Crossing over
    if params.crossing_over
        # Shuffle the order of the Genomes.
        shuffled_genomes = shuffle(child_genomes)
        child_genomes = []
        # Perform crossing over by pairs; ignore the last Genome if
        # there's an odd number.
        for i in 1:div(length(shuffled_genomes), 2)
            genome1, genome2 = cross_over(shuffled_genomes[i*2-1],
                                          shuffled_genomes[i*2])
            push!(child_genomes, genome1)
            push!(child_genomes, genome2)
        end
    end

    return child_genomes
end
