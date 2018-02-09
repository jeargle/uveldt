# John Eargle (mailto: jeargle at gmail.com)
# 2018
# Uveldt

module Uveldt

export Element, ElementTable, Molecule, add_element


type Element
    name::AbstractString
    mass::Int64
    Element(name::AbstractString, mass::Int64) = new(name, mass)
end

# ElementTable acts like a periodic table.
type ElementTable
    elements::Dict{AbstractString, Element}
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

# Add an Element to an ElementTable.
# element_table: ElementTable that receives the Element
# element: Element to add
function add_element(element_table::ElementTable, element::Element)
    element_table.elements[element.name] = element
end


end
