# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt.evolution

"""
Substitution matrix for SNPs.
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
    snp_rate::Float64
    insertion_rate::Float64
    deletion_rate::Float64
    duplication_rate::Float64
    inversion_rate::Float64
    translocation_rate::Float64
    substitution_matrix::Union{SubstitutionMatrix, Nothing}

    MutationParams(snp_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate) = new(snp_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, nothing)
    MutationParams(snp_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, substitution_matrix) = new(snp_rate, insertion_rate, deletion_rate, duplication_rate, inversion_rate, translocation_rate, substitution_matrix)
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
column j.  SNP target selection already takes into account which bases
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
    add_snps(genome, rate; sub_mat)

Add SNPs to Genome.

# Arguments
- genome::Genome
- rate
- sub_mat::SubstitutionMatrix

# Returns
- Mutated Genome
"""
function add_snps(genome::Genome, rate; sub_mat::Union{SubstitutionMatrix, Nothing}=nothing)
    geom_dist = Geometric(rate)
    location = 1 + rand(geom_dist)
    alphabet = alphabet_string(genome.chemistry)
    fragments = []
    start = 1

    while location <= length(genome.string)
        push!(fragments, genome.string[start:location-1])
        if sub_mat == nothing
            snp = string(rand(alphabet))
        else
            snp = string(substitute(sub_mat, genome.string[location]))
        end
        @printf "  %d SNP %s->%s\n" location genome.string[location] snp
        push!(fragments, snp)
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
        genome = add_snps(genome, params.snp_rate; sub_mat=params.substitution_matrix)
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
