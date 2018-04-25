# John Eargle (mailto: jeargle at gmail.com)
# 2018
# uveldt

module uveldt

export Element, ElementTable, Bond, BondTable, Chemistry, Molecule, Reaction, Gene, Genome, Cell, VeldtPoint, Veldt, init_molecules, parse_reaction, transcribe_gene, genome_string, find_genes, is_pseudogene, read_fasta, write_fasta, get_bond, mass


type Element
    name::Char
    mass::Int64
    Element(name::Char, mass::Int64) = new(name, mass)
end

Base.show(io::IO, el::Element) = show(io, string(el.name))
Base.show(io::IO, m::MIME"text/plain", el::Element) = show(io, m, string(el.name))


"""
ElementTable acts like a periodic table.
"""
type ElementTable
    elements::Dict{Char, Element}

    function ElementTable(elements)
        element_dict = Dict()
        for element in elements
            element_dict[element.name] = element
        end
        new(element_dict)
    end
end

Base.show(io::IO, et::ElementTable) = show(io, keys(et.elements))
Base.show(io::IO, m::MIME"text/plain", et::ElementTable) = show(io, m, keys(et.elements))


"""
Bond between two Elements
"""
type Bond
    element1::Element
    element2::Element
    element_table::ElementTable
    energy_change::Float64
    transition_energy::Float64

    Bond(element1::Element,
         element2::Element,
         element_table::ElementTable,
         energy_change::Float64,
         transition_energy::Float64) = new(element1, element2, element_table,
                                           energy_change, transition_energy)
end

Base.show(io::IO, b::Bond) = show(io, string(b.element1.name, '-', b.element2.name, ", trans: ", b.transition_energy, ", energy: ", b.energy_change))
Base.show(io::IO, m::MIME"text/plain", b::Bond) = show(io, m, string(b.element1.name, "-", b.element2.name, ", trans: ", b.transition_energy, ", energy: ", b.energy_change))


"""
Set of all Bonds
"""
type BondTable
    bonds::Dict{Char, Dict{Char, Bond}}
    element_table::ElementTable

    function BondTable(bonds, element_table::ElementTable)
        bond_dict = Dict{Char, Dict{Char, Bond}}()

        for bond in bonds
            if haskey(bond_dict, bond.element1.name)
                bond_dict[bond.element1.name][bond.element2.name] = bond
            else
                bond_dict[bond.element1.name] = Dict(bond.element2.name=>bond)
            end

            # Make sure entries are added for both lookup orders
            if bond.element1 != bond.element2
                if haskey(bond_dict, bond.element2.name)
                    bond_dict[bond.element2.name][bond.element1.name] = bond
                else
                    bond_dict[bond.element2.name] = Dict(bond.element1.name=>bond)
                end
            end
        end

        new(bond_dict, element_table)
    end
end

Base.show(io::IO, bt::BondTable) = show(io, string(bt.bonds))
Base.show(io::IO, m::MIME"text/plain", bt::BondTable) = show(io, m, string(bt.bonds))


"""
Container for Elements and Bonds.
"""
type Chemistry
    element_table::ElementTable
    bond_table::BondTable

    function Chemistry(element_table::ElementTable, bond_table::BondTable)
        # check bond_table for full set of Bonds for all Element pairs
        new(element_table, bond_table)
    end
end


"""
1D String of Elements.
"""
type Molecule
    name::AbstractString
    elements::AbstractString
    chemistry::Chemistry

    Molecule(name::AbstractString,
             elements::AbstractString,
             chemistry::Chemistry) = new(name, elements, chemistry)
end

Base.show(io::IO, mol::Molecule) = show(io, string(mol.name, ", ", mol.elements))
Base.show(io::IO, m::MIME"text/plain", mol::Molecule) = show(io, m, string(mol.name, ", ", mol.elements))


"""
Molecular reaction with reactants and products where Bonds are created and/or broken.  A
Reaction is specified by a string derived from a Gene.
"""
type Reaction
    reactants::Array{Molecule, 1}
    products::Array{Molecule, 1}
    old_bonds::Array{Bond, 1}
    new_bonds::Array{Bond, 1}
    energy_change::Float64
    transition_energy::Float64

    function Reaction(reactants, products, chemistry::Chemistry)
        energy_change = 0.0
        transition_energy = 0.0

        r_mols = []
        for r_str in reactants
            mol = Molecule(r_str, r_str, chemistry)
            push!(r_mols, mol)
        end

        new_bonds = []
        for i in 1:length(reactants)-1
            el1 = reactants[i][end]
            el2 = reactants[i+1][1]
            bond = get_bond(chemistry.bond_table, el1, el2)
            push!(new_bonds, bond)
            energy_change += bond.energy_change
            transition_energy += bond.transition_energy
        end

        p_mols = []
        for p_str in products
            mol = Molecule(p_str, p_str, chemistry)
            push!(p_mols, mol)
        end

        old_bonds = []
        for i in 1:length(products)-1
            el1 = products[i][end]
            el2 = products[i+1][1]
            bond = get_bond(chemistry.bond_table, el1, el2)
            push!(old_bonds, bond)
            energy_change -= bond.energy_change
            transition_energy += bond.transition_energy - bond.energy_change
        end

        new(r_mols, p_mols, old_bonds, new_bonds,
            energy_change, transition_energy)
        # new(r_mols, p_mols, old_bonds, new_bonds, energy_change, transition_energy)
    end
end

Base.show(io::IO, r::Reaction) = show(io, "reactants: " * join([m.name for m in r.reactants], ", ") * " products: " * join([m.name for m in r.products], ", ") * ", trans: " * string(r.transition_energy) * ", energy: " * string(r.energy_change))
Base.show(io::IO, m::MIME"text/plain", r::Reaction) = show(io, m, "reactants: " * join([m.name for m in r.reactants], ", ") * " products: " * join([m.name for m in r.products], ", ") * ", trans: " * string(r.transition_energy) * ", energy: " * string(r.energy_change))


"""
Genes are strings bracketed by start '(' and stop ')' characters that code for specific
chemical reactions.  The internals consist of Molecule strings with join '*' and split '/'
operators specifying how the set of Molecules will be processed.  All operators act at the
same time.  Each join creates a new bond, and each split breaks an existing bond.

# Examples
  (A*B): A + B -> AB
  (A/B): AB -> A + B
  (A*B/C): A + BC -> AB + C
  (A*B*C): A + B + C -> ABC
  (A*A*A): A + A + A -> AAA
"""
type Gene
    location::Int64
    string::AbstractString
    transcript::AbstractString
    element_table::ElementTable
    function Gene(location::Int64, gene_string::AbstractString, element_table::ElementTable)
        transcript = transcribe_gene(gene_string, element_table)
        new(location, gene_string, transcript, element_table)
    end
end


"""
String specifying all genes in the cell.  Genes are regions of the Genome that fit the Gene
pattern.  A Genome can have a lot of intergenic space that doesn't code for anything.
"""
type Genome
    name::AbstractString
    string::AbstractString
    element_table::ElementTable
    genes::Array{Gene, 1}

    Genome(name::AbstractString, string::AbstractString, element_table::ElementTable) = new(name, string, element_table)

    function Genome(name::AbstractString, size::Int64, element_table::ElementTable)
        Genome(name, genome_string(size, element_table), element_table)
    end
end


"""
Membrane-bound compartment containing Molecules and a Genome.  Chemical Reactions happen
within.
"""
type Cell
    genome::Genome
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    energy::Int64

    function Cell(genome::Genome)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(2)
        molecule_counts[1] = Dict{AbstractString, Int64}()
        molecule_counts[2] = Dict{AbstractString, Int64}()
        new(genome, molecule_counts, 0)
    end
end


"""
Unique location in a Veldt.  Contains two Molecule count buffers and possibly a single Cell.
VeldtPoint time evolution is controlled by a Veldt.
"""
type VeldtPoint
    molecule_counts::Array{Dict{AbstractString, Int64}, 1}  # 2 Dicts: current, future
    cell::Cell

    # function VeldtPoint()
    #     molecule_counts = Array{Dict{AbstractString, Int64}, 1}(2)
    #     molecule_counts[1] = Dict{AbstractString, Int64}()
    #     molecule_counts[2] = Dict{AbstractString, Int64}()
    #     new(molecule_counts)
    # end

    function VeldtPoint(; molecules=[], cell=nothing)
        molecule_counts = Array{Dict{AbstractString, Int64}, 1}(2)
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
type Veldt
    dims::Array{Int64, 1}   # 2 or 3
    points  # Array with dims dimensions holding VeldtPoints
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
        new(dims, points, Dict{AbstractString, Int64}(), Dict{AbstractString, Int64}())
    end
end


"""
    init_molecules()

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
    parse_reaction(reaction_string, element_table)

Parse the reactants and products from a reaction string.
"""
function parse_reaction(reaction_string::AbstractString, element_table::ElementTable)
    reactants = split(replace(reaction_string, "/", ""), "*")
    products = split(replace(reaction_string, "*", ""), "/")
    return (reactants, products)
end

"""
    transcribe_gene(gene_string, element_table)

Convert a Gene string into a valid Reaction string by removing redundant directive
characters.
"""
function transcribe_gene(gene_string::AbstractString, element_table::ElementTable)
    # modes for parsing
    #   0: in Molecule
    #   1: in directive
    mode = 0
    gene_chars = []
    elements = join(keys(element_table.elements))

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
            elseif search(elements, el) > 0
                push!(gene_chars, el)
                mode = 0
            end
        elseif mode == 1
            if search(elements, el) > 0
                push!(gene_chars, el)
                mode = 0
            elseif el != '*' && el != '/'
                @printf("Error: bad character in Gene string - %s", el)
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
    genome_string(size, element_table)

Randomly generate a valid Genome string.
"""
function genome_string(size::Int64, element_table::ElementTable)
    genome = Array{AbstractString}(size)
    elements = [string(el) for el in keys(element_table.elements)]
    rand!(genome, vcat(["(", ")", "*", "/"], elements))
    return join(genome)
end

"""
    find_genes(genome)

Search through a Genome for patterns that match possible Gene strings and build Genes from
them.
"""
function find_genes(genome::Genome)
    gene_match = eachmatch(r"\([^\(]*?\)", genome.string)
    genes = [Gene(gm.offset, gm.match, genome.element_table)
             for gm in collect(gene_match)]
    return genes
end


"""
    is_pseudogene(gene)

Identify pseudogene.
"""
function is_pseudogene(gene::Gene)
    elements = join(keys(gene.element_table.elements))
    if length(gene.transcript) < 3
        return true
    end

    m = match(Regex("([$(elements)].*[$(elements)])"), gene.transcript)
    if m === nothing
        return true
    end

    return !ismatch(r"[\*/]", m[1])

    # return !ismatch(Regex("[$(elements)]+"), gene.transcript)
end


"""
    read_fasta(filename)

Read a FASTA file and return the information as tuples of names and a genome strings.
"""
function read_fasta(filename)
    genome_info = []
    name = ""
    genome_str = ""
    first_genome = true
    open(filename, "r") do f
        for line in eachline(f)
            print(line)
            if line[1] == '>'
                if !first_genome
                    push!(genome_info, (name, genome_str))
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
    write_fasta(genomes)

Write a multi-FASTA file with the full genome strings for an array of Genomes.
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
    write_fasta(genome)

Write a FASTA file with the full genome string.
"""
function write_fasta(genome::Genome, filename)
    write_fasta([genome], filename)
end


"""
    get_bond(bond_table, element1, element2)

Add a Bond to the a BondTable.
"""
function get_bond(bond_table::BondTable, element1::Char, element2::Char)
    return bond_table.bonds[element1][element2]
end


"""
    mass(element)

Get the mass of an Element.
"""
function mass(element::Element)
    return element.mass
end


"""
    mass(molecule)

Get the mass of a Molecule.
"""
function mass(molecule::Molecule)
    total_mass = 0
    for el in molecule.elements
        total_mass += molecule.chemistry.element_table.elements[el].mass
    end
    return total_mass
end


end
