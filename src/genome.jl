# John Eargle (mailto: jeargle at gmail.com)
# uveldt.genome

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
    uuid::UUID
    name::AbstractString
    string::AbstractString
    chemistry::Chemistry
    genes::Array{Gene, 1}

    Genome(uuid::UUID, name::AbstractString, string::AbstractString, chemistry::Chemistry) = new(uuid, name, string, chemistry)

    Genome(name::AbstractString, string::AbstractString, chemistry::Chemistry) = new(UUIDs.uuid4(), name, string, chemistry)

    Genome(string::AbstractString, chemistry::Chemistry) = new(UUIDs.uuid4(), "", string, chemistry)

    function Genome(uuid::UUID, name::AbstractString, size::Int64, chemistry::Chemistry)
        Genome(uuid, name, genome_string(size, chemistry), chemistry)
    end

    function Genome(name::AbstractString, size::Int64, chemistry::Chemistry)
        Genome(name, genome_string(size, chemistry), chemistry)
    end
end


"""
    parse_reaction(reaction_string::AbstractString)

Parse the reactants and products from a reaction string.

# Arguments
- `reaction_string::AbstractString`:

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
    transcribe_gene(gene_string::AbstractString, chemistry::Chemistry)

Convert a Gene string into a valid Reaction string by removing redundant directive
characters.

# Arguments
- `gene_string::AbstractString`:
- `chemistry::Chemistry`:

# Returns
- `AbstractString`: Reaction string
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

    # println("  elements: ", elements)
    # println("  in: ", gene_string)

    for el in gene_string[2:end-1]
        if mode == 0
            if el == '*' || el == '/'
                push!(gene_chars, el)
                mode = 1
            elseif in(el, elements)
                push!(gene_chars, el)
                mode = 0
            end
        elseif mode == 1
            if in(el, elements)
                push!(gene_chars, el)
                mode = 0
            elseif el != '*' && el != '/'
                @printf "Error: bad character in Gene string - %s\n" el
            end
        end
    end

    # Do not end string with a directive character.
    if mode == 1
        pop!(gene_chars)
    end

    # println("  out: ", join(gene_chars))
    return join(gene_chars)
end


"""
    genome_string(size::Int64, chemistry::Chemistry)

Randomly generate a valid Genome string.

# Arguments
- `size::Int64`:
- `chemistry::Chemistry`:

# Returns
- `AbstractString`: Genome string
"""
function genome_string(size::Int64, chemistry::Chemistry)
    return randstring(alphabet_string(chemistry), size)
end


"""
    find_genes(genome::Genome)

Search through a Genome for patterns that match possible Gene strings and build Genes from
them.

# Arguments
- `genome::Genome`:

# Returns
- ``: Array of Genes
"""
function find_genes(genome::Genome)
    rx = r"\([^\(]*?\)"
    gene_match = eachmatch(rx, genome.string)
    genes = [Gene(gm.offset, gm.match, genome.chemistry)
             for gm in gene_match]
    return genes
end


"""
    is_pseudogene(gene::Gene)

Identify pseudogene.

There currently aren't many ways to code a pseudogene.  The main one
is to have an empty transcript.  Repeated and tail-end directives get
ignored.

# Arguments
- `gene::Gene`:

# Returns
- `Bool`:
"""
function is_pseudogene(gene::Gene)
    elements = join(keys(gene.chemistry.element_table.elements))

    if length(gene.transcript) == 0
        return true
    end

    m = match(r"([\*/]*[$(elements)]+[\*/$(elements)]*)", gene.string)
    if m === nothing
        return true
    end

    return false
end


"""
    read_fasta(filename)

Read a FASTA file and return the information as tuples of names and a genome strings.

# Arguments
- `filename`: name of input FASTA file

# Returns
- ``: Array of tuples (name string, genome string)
"""
function read_fasta(filename)
    genome_info = []
    name = ""
    uuid = nothing
    genome_str = ""
    first_genome = true
    open(filename, "r") do f
        for line in eachline(f)
            println(line)
            if line[1] == '>'
                if !first_genome
                    push!(genome_info, (name, uuid, genome_str))
                else
                    first_genome = false
                end
                id_info = split(line[2:end], " ")
                name = id_info[1]
                if length(id_info) == 1
                    uuid = nothing
                else
                    uuid = UUID(id_info[2])
                end
                genome_str = ""
            else
                genome_str *= line
            end
        end
    end
    push!(genome_info, (name, uuid, genome_str))
    return genome_info
end


"""
    write_fasta(genomes, filename)

Write a multi-FASTA file with the full genome strings for an array of Genomes.

# Arguments
- `genomes`: array of Genomes
- `filename`: name of output FASTA file
"""
function write_fasta(genomes, filename)
    f = open(filename, "w")
    line_length = 80
    for genome in genomes
        println(f, ">", genome.name, " ", genome.uuid)
        full_line_count = div(length(genome.string), line_length)  # integer division
        for i in 1:full_line_count
            println(f, genome.string[(i-1)*line_length+1:i*line_length])
        end
        println(f, genome.string[full_line_count*line_length+1:length(genome.string)])
    end
    close(f)
end


"""
    write_fasta(genome::Genome, filename)

Write a FASTA file with the genome string for a single Genome.

# Arguments
- `genome::Genome`: Genome
- `filename`: name of output FASTA file
"""
function write_fasta(genome::Genome, filename)
    write_fasta([genome], filename)
end
