# John Eargle (mailto: jeargle at gmail.com)
# uveldt.genome

@exported_enum GenomeStrand positive negative

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
    strand::GenomeStrand

    function Gene(location::Int64, gene_string::AbstractString, chemistry::Chemistry, strand::GenomeStrand)
        transcript = transcribe_gene(gene_string, chemistry)
        new(location, gene_string, transcript, chemistry, strand)
    end

    function Gene(location::Int64, gene_string::AbstractString, chemistry::Chemistry)
        transcript = transcribe_gene(gene_string, chemistry)
        new(location, gene_string, transcript, chemistry, positive)
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
    parse_reaction(reaction_string::AbstractString, chemistry::Chemistry)

Parse the reactants and products from a reaction string.

# Arguments
- `reaction_string::AbstractString`: string representation of a Reaction
- `chemistry::Chemistry`: Chemistry for Reactions to use

# Returns
- `Reaction`
"""
function parse_reaction(reaction_string::AbstractString, chemistry::Chemistry)
    if isnothing(match(r"[\*/]", reaction_string))
        reaction_type = pore
        reaction1 = Reaction(reaction_string, chemistry)
    elseif reaction_string[1] == '*'
        reaction_type = transport_in
        reactant = replace(reaction_string, r"\*|/" => "")
        reaction1 = Reaction(reactant, chemistry, reaction_type, -3.0)
    elseif reaction_string[1] == '/'
        reaction_type = transport_out
        reactant = replace(reaction_string, r"\*|/" => "")
        reaction1 = Reaction(reactant, chemistry, reaction_type, -3.0)
    else
        reactants = split(replace(reaction_string, "/" => ""), "*")
        products = split(replace(reaction_string, "*" => ""), "/")
        reaction1 = Reaction(reactants, products, chemistry)
    end

    return reaction1
end


"""
    parse_reactions(genome::Genome)

Search through a Genome and find all of the Reactions it codes for.

# Arguments
- `genome::Genome`: Genome from which to extract Reactions

# Returns
- `Array{Reaction, 1}`: Array of Reactions
"""
function parse_reactions(genome::Genome)
    println("  * Gene transcripts")
    genes = find_genes(genome)
    reactions = Array{Reaction, 1}()
    chemistry = genome.chemistry

    for gene in genes
        # @printf "gene%d.transcript: %s\n" i genes[i].transcript
        if !is_pseudogene(gene)
            push!(reactions, parse_reaction(gene.transcript))
        end
    end

    return reactions
end


"""
    transcribe_gene(gene_string::AbstractString, chemistry::Chemistry)

Convert a Gene string into a valid Reaction string by removing redundant directive
characters.

# Arguments
- `gene_string::AbstractString`: string representation of a Gene
- `chemistry::Chemistry`: Chemistry for Reactions to use

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
- `size::Int64`: size, in characters, that the Genome should be
- `chemistry::Chemistry`: Chemistry for Genome to use

# Returns
- `AbstractString`: string representation of Genome contents
"""
function genome_string(size::Int64, chemistry::Chemistry)
    return randstring(alphabet_string(chemistry), size)
end


"""
    find_genes(genome_string::AbstractString, chemistry::Chemistry)

Search through a genome string for patterns that match possible Gene strings and build
Genes from them.

# Arguments
- `genome_string::AbstractString`: string representation of Genome contents
- `chemistry::Chemistry`:  Chemistry for Genes to use
- `strand::GenomeStrand=positive`: positive of negative

# Returns
- `Array{Gene, 1}`: Array of Genes
"""
function find_genes(genome_string::AbstractString, chemistry::Chemistry; strand::GenomeStrand=positive)
    rx = r"\([^\(]*?\)"
    gene_match = eachmatch(rx, genome_string)
    genes = [Gene(gm.offset, gm.match, chemistry, strand)
             for gm in gene_match]
    return genes
end


"""
    find_genes(genome::Genome; include_reverse::Bool=false)

Search through a Genome for patterns that match possible Gene strings and build Genes from
them.

# Arguments
- `genome::Genome`:  Genome from which to extract Genes
- `include_reverse::Bool=false`: whether or not to find genes by reading the Genome backwards

# Returns
- `Array{Gene, 1}`: Array of Genes
"""
function find_genes(genome::Genome; include_reverse::Bool=false)
    genes = find_genes(genome.string, genome.chemistry)

    if include_reverse
        genes_reverse = find_genes(reverse(genome.string), genome.chemistry, strand=negative)
        genes = vcat(genes, genes_reverse)
    end

    return genes
end


"""
    is_pseudogene(gene::Gene)

Identify pseudogene.

There currently aren't many ways to code a pseudogene.  The main one
is to have an empty transcript.  Repeated and tail-end directives get
ignored.

# Arguments
- `gene::Gene`: Gene to examine

# Returns
- `Bool`: whether or not the Gene is a pseudogene
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
- `Array{(AbstractString, AbstractString), 1}`: Array of tuples (name string, genome string)
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
- `genome::Genome`: Genome to write to a FASTA file
- `filename`: name of output FASTA file
"""
function write_fasta(genome::Genome, filename)
    write_fasta([genome], filename)
end
