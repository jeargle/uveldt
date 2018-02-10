# John Eargle (mailto: jeargle at gmail.com)
# 2018
# Uveldt

module Uveldt

export Element, ElementTable, Molecule, Bond, BondTable, add_element, add_bond, get_bond, mass


type Element
    name::Char
    mass::Int64
    Element(name::Char, mass::Int64) = new(name, mass)
end

# ElementTable acts like a periodic table.
type ElementTable
    elements::Dict{Char, Element}
    ElementTable() = new(Dict())
end

# 1D String of Elements.
type Molecule
    name::AbstractString
    elements::AbstractString
    element_table::ElementTable
    Molecule(name::AbstractString,
             elements::AbstractString,
             element_table::ElementTable) = new(name, elements, element_table)
end

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

# Set of all Bonds
type BondTable
    bonds::Dict{Char, Dict{Char, Bond}}
    element_table::ElementTable
    BondTable(element_table::ElementTable) = new(Dict(), element_table)
end


# Add an Element to an ElementTable.
# element_table: ElementTable that receives the Element
# element: Element to add
function add_element(element_table::ElementTable, element::Element)
    element_table.elements[element.name] = element
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
