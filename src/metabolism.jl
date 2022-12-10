# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt.metabolism


"""
Graph with nodes representing Molecules and Reactions.
"""
struct Pathway
    molecules::Array{Molecule, 1}
    reactions::Array{Reaction, 1}
    reactant_edges::Array
    product_edges::Array

    function Pathway(molecules::Array{Molecule, 1},
                     reactions::Array{Reaction, 1},
                     reactant_edges::Array,
                     product_edges::Array)
        new(molecules, reactions, reactant_edges, product_edges)
    end
end


"""
Set of disconnected Pathways.
"""
struct Metabolism
    pathways::Array{Pathway, 1}

    function Metabolism(pathways::Array{Pathway, 1})
        new(pathways)
    end
end
