# John Eargle (mailto: jeargle at gmail.com)
# 2018
# Uveldt

module Uveldt

export Element, ElementTable, Molecule, Bond, add_element, mass


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


# Add an Element to an ElementTable.
# element_table: ElementTable that receives the Element
# element: Element to add
function add_element(element_table::ElementTable, element::Element)
    element_table.elements[element.name] = element
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
