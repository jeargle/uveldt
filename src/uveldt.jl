# John Eargle (mailto: jeargle at gmail.com)
# 2018
# Uveldt

module Uveldt

export Element, ElementTable, Molecule, Bond, BondTable, Gene, Genome, build_genome, find_genes, add_elements, add_bond, get_bond, mass


type Element
    name::Char
    mass::Int64
    Element(name::Char, mass::Int64) = new(name, mass)
end

Base.show(io::IO, el::Element) = show(io, string(el.name))
Base.show(io::IO, m::MIME"text/plain", el::Element) = show(io, m, string(el.name))


# ElementTable acts like a periodic table.
type ElementTable
    elements::Dict{Char, Element}
    ElementTable() = new(Dict())
end

Base.show(io::IO, et::ElementTable) = show(io, keys(et.elements))
Base.show(io::IO, m::MIME"text/plain", et::ElementTable) = show(io, m, keys(et.elements))


# 1D String of Elements.
type Molecule
    name::AbstractString
    elements::AbstractString
    element_table::ElementTable
    Molecule(name::AbstractString,
             elements::AbstractString,
             element_table::ElementTable) = new(name, elements, element_table)
end

Base.show(io::IO, mol::Molecule) = show(io, string(mol.name, ", ", mol.elements))
Base.show(io::IO, m::MIME"text/plain", mol::Molecule) = show(io, m, string(mol.name, ", ", mol.elements))


# Bond between two Elements
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


# Set of all Bonds
type BondTable
    bonds::Dict{Char, Dict{Char, Bond}}
    element_table::ElementTable
    BondTable(element_table::ElementTable) = new(Dict(), element_table)
end

Base.show(io::IO, bt::BondTable) = show(io, string(bt.bonds))
Base.show(io::IO, m::MIME"text/plain", bt::BondTable) = show(io, m, string(bt.bonds))


# Genes are strings bracketed by start '(' and stop ')' characters
# that code for specific chemical reactions.  The internals consist
# of Molecule strings with join '*' and split '/' operators specifying
# how the set of Molecules will be processed.  All operators act at
# the same time.  Each join creates a new bond, and each split breaks
# an existing bond.
#
# Examples
#   (A*B): A + B -> AB
#   (A/B): AB -> A + B
#   (A*B/C): A + BC -> AB + C
#   (A*B*C): A + B + C -> ABC
#   (A*A*A): A + A + A -> AAA
type Gene
    location::Int64
    gene::AbstractString
    Gene(location::Int64, gene::AbstractString) = new(location, gene)
end


# String specifying all genes in the cell.  Genes are regions of the
# Genome that fit the Gene pattern.  A Genome can have a lot of
# intergenic space that doesn't code for anything.
type Genome
    name::AbstractString
    genome::AbstractString
    element_table::ElementTable
    genes::Array{Gene, 1}
    Genome(name::AbstractString, genome::AbstractString, element_table::ElementTable) = new(name, genome, element_table)
end

function build_genome(size::Int64, element_table::ElementTable)
    genome = Array{AbstractString}(size)
    elements = [string(el) for el in keys(element_table.elements)]
    rand!(genome, vcat(["(", ")", "*", "/"], elements))
    return join(genome)
end

# Search through a Genome for patterns that match possible Gene
# strings and build Genes from them.
function find_genes(genome::Genome)
    gene_match = eachmatch(r"\([^\(]*?\)", genome.genome)
    # genes = [gm.match for gm in gene_match]
    for gene in gene_match
        println(gene)
    end
end

# Write a FASTA file with the full genome string.
function write_fasta(genome::Genome)
end


# Add an Element to an ElementTable.
# element_table: ElementTable that receives the Element
# element: Element to add
function add_elements(element_table::ElementTable, elements::Array{Element, 1})
    for element in elements
        element_table.elements[element.name] = element
    end
end

# Add a Bond to the a BondTable.
# bond_table: BondTable that receives the Bond
# bond: Bond to add
function add_bond(bond_table::BondTable, bond::Bond)
    if haskey(bond_table.bonds, bond.element1.name)
        bond_table.bonds[bond.element1.name][bond.element2.name] = bond
    else
        bond_table.bonds[bond.element1.name] = Dict(bond.element2.name=>bond)
    end

    # Make sure entries are added for both lookup orders
    if bond.element1 != bond.element2
        if haskey(bond_table.bonds, bond.element2.name)
            bond_table.bonds[bond.element2.name][bond.element1.name] = bond
        else
            bond_table.bonds[bond.element2.name] = Dict(bond.element1.name=>bond)
        end
    end
end

# Add a Bond to the a BondTable.
# bond_table: BondTable
# element1: name of first Element
# element2: name of second Element
function get_bond(bond_table::BondTable, element1::Char, element2::Char)
    return bond_table.bonds[element1][element2]
end

function mass(element::Element)
    return element.mass
end

function mass(molecule::Molecule)
    total_mass = 0
    for el in molecule.elements
        total_mass += molecule.element_table.elements[el].mass
    end
    return total_mass
end


end
