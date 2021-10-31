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
Evolution parameters for evolution algorithm.
"""
struct MutationParams
    snv_rate::Float64
    insertion_rate::Float64
    deletion_rate::Float64
    duplication_rate::Float64
    inversion_rate::Float64
    translocation_rate::Float64
    substitution_matrix::Union{SubstitutionMatrix, Nothing}

    MutationParams(snv_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate) = new(snv_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, nothing)
    MutationParams(snv_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, substitution_matrix) = new(snv_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, substitution_matrix)
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
        @printf "  %d SNV %s->%s\n" location genome.string[location] snv
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
        @printf "  %d insert %s\n" location insert
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
        @printf "  %d delete %s\n" location genome.string[location:location+size]
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
        @printf "  %d duplicate %s\n" location duplication
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
    @printf "  length(genome.string) %d\n" length(genome.string)
    for (insertion_point, duplication) in ordered_inserts
        @printf "  start: %d insertion_point-1 %d duplication %s\n" start insertion_point-1 duplication
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
    mutate(genomes, params)

Take an Array of parent Genomes and create a set of child Genomes.

# Arguments
- genomes::Array{Genome, 1}

# Returns
- Array{Genome, 1} containing new child Genomes
"""
function mutate(genomes, params::MutationParams)
    child_genomes = []

    for genome in genomes
        genome = add_snvs(genome, params.snv_rate; sub_mat=params.substitution_matrix)
        genome = add_insertions(genome, params.insertion_rate; size_param=0.5)
        genome = remove_deletions(genome, params.deletion_rate; size_param=0.5)
        # # Duplications
        # add_duplications(genome, params.duplication_rate; size_param=0.5)
        # # Inversions
        # add_inversions(genome, params.inversion_rate; size_param=0.5)
        # # Translocations
        # add_translocations(genome, params.translocation_rate; size_param=0.5)
        push!(child_genomes, genome)
    end

    # # Crossing over
    # cross_over(genome1, genome2)

    return child_genomes
end
