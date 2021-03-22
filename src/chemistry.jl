# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt

macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end


"""
Chemical element.
"""
struct Element
    name::Char
    mass::Int64
    Element(name::Char, mass::Int64) = new(name, mass)
end

Base.show(io::IO, el::Element) = show(io, string(el.name))
Base.show(io::IO, m::MIME"text/plain", el::Element) = show(io, m, string(el.name))


"""
Periodic table for Elements.
"""
struct ElementTable
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
Bond between two Elements.
"""
struct Bond
    element1::Element
    element2::Element
    energy_change::Float64
    transition_energy::Float64

    Bond(element1::Element,
         element2::Element,
         energy_change::Float64,
         transition_energy::Float64) = new(element1, element2, energy_change,
                                           transition_energy)
end

Base.show(io::IO, b::Bond) = show(io, string(b.element1.name, '-', b.element2.name, ", trans: ", b.transition_energy, ", energy: ", b.energy_change))
Base.show(io::IO, m::MIME"text/plain", b::Bond) = show(io, m, string(b.element1.name, "-", b.element2.name, ", trans: ", b.transition_energy, ", energy: ", b.energy_change))


"""
Set of all Bonds.
"""
struct BondTable
    bonds::Dict{Char, Dict{Char, Bond}}

    function BondTable(bonds)
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

        new(bond_dict)
    end
end

Base.show(io::IO, bt::BondTable) = show(io, string(bt.bonds))
Base.show(io::IO, m::MIME"text/plain", bt::BondTable) = show(io, m, string(bt.bonds))


"""
Container for Elements and Bonds.
"""
struct Chemistry
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
struct Molecule
    name::AbstractString
    elements::AbstractString
    chemistry::Chemistry

    Molecule(name::AbstractString,
             elements::AbstractString,
             chemistry::Chemistry) = new(name, elements, chemistry)
end

Base.show(io::IO, mol::Molecule) = show(io, string(mol.name, ", ", mol.elements))
Base.show(io::IO, m::MIME"text/plain", mol::Molecule) = show(io, m, string(mol.name, ", ", mol.elements))


@exported_enum ReactionType reaction pore transport_in transport_out

"""
Molecular reaction with reactants and products where Bonds are created and/or broken.  A
Reaction is specified by a string derived from a Gene.
"""
struct Reaction
    reaction_type::ReactionType
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

        new(reaction, r_mols, p_mols, old_bonds, new_bonds,
            energy_change, transition_energy)
        # new(r_mols, p_mols, old_bonds, new_bonds, energy_change, transition_energy)
    end

    function Reaction(reactant, chemistry::Chemistry)
        energy_change = 0.0
        transition_energy = 0.0

        r_mol = Molecule(reactant, reactant, chemistry)

        new(pore, [r_mol], [], [], [],
            energy_change, transition_energy)
    end

    function Reaction(reactant, chemistry::Chemistry,
                      transport_type::ReactionType, energy_change)
        transition_energy = 0.0

        if !(transport_type in (transport_in, transport_out))
            println("Error: bad transport_type")
        end

        r_mol = Molecule(reactant, reactant, chemistry)

        new(transport_type, [r_mol], [], [], [],
            energy_change, transition_energy)
    end
end

Base.show(io::IO, r::Reaction) = show(io, string(r.reaction_type) * " reactants: " * join([m.name for m in r.reactants], ", ") * " products: " * join([m.name for m in r.products], ", ") * ", trans: " * string(r.transition_energy) * ", energy: " * string(r.energy_change))
Base.show(io::IO, m::MIME"text/plain", r::Reaction) = show(io, m, string(r.reaction_type) * "reactants: " * join([m.name for m in r.reactants], ", ") * " products: " * join([m.name for m in r.products], ", ") * ", trans: " * string(r.transition_energy) * ", energy: " * string(r.energy_change))


"""
    alphabet_string(chemistry::Chemistry)

Generate a string with all possible characters that can appear with
this Chemistry.

# Arguments
- chemistry::Chemistry

# Returns
- String
"""
function alphabet_string(chemistry::Chemistry)
    elements = [el for el in keys(chemistry.element_table.elements)]
    return String(elements) * "()*/"
end


"""
    get_bond(bond_table, element1, element2)

Get a Bond from a BondTable.

# Arguments
- bond_table::BondTable
- element1::Char
- element2::Char

# Returns
- Bond
"""
function get_bond(bond_table::BondTable, element1::Char, element2::Char)
    return bond_table.bonds[element1][element2]
end


"""
    mass(element)

Get the mass of an Element.

# Arguments
- element::Element

# Returns
- Mass of the Elment
"""
function mass(element::Element)
    return element.mass
end


"""
    mass(molecule::Molecule)

Get the mass of a Molecule.

# Arguments
- molecule::Molecule

# Returns
- Mass of the Molecule
"""
function mass(molecule::Molecule)
    total_mass = 0
    for el in molecule.elements
        total_mass += molecule.chemistry.element_table.elements[el].mass
    end
    return total_mass
end


"""
    setup_chemistry(filename)

Create a Chemistry from a YAML setup file.

# Arguments
- filename: name of YAML setup file

# Returns
- Chemistry
"""
function setup_chemistry(filename)
    setup = YAML.load(open(filename))

    # build Elements
    elements = Array{Element, 1}()

    if haskey(setup, "elements")
        for el_info in setup["elements"]
            name = el_info["name"][1]
            mass = el_info["mass"]
            element = Element(name, mass)
            push!(elements, element)
        end
    end

    element_table = ElementTable(elements)

    # build Bonds
    bonds = Array{Bond, 1}()

    if haskey(setup, "bonds")
        for bond_info in setup["bonds"]
            el1 = element_table.elements[bond_info["elements"][1][1]]
            el2 = element_table.elements[bond_info["elements"][2][1]]
            energy_change = bond_info["energy_change"]
            transition_energy = bond_info["transition_energy"]
            bond = Bond(el1, el2, energy_change, transition_energy)
            push!(bonds, bond)
        end
    end

    bond_table = BondTable(bonds)

    return (Chemistry(element_table, bond_table))
end
