# John Eargle (mailto: jeargle at gmail.com)
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

If nodes or edges are present in multiple Pathways, they will be merged
into a connected graph.  The individual Pathway representations will be
maintained, but the base graph will be the result of all possible merges
between Pathways.
"""
struct Metabolism
    pathways::Array{Pathway, 1}
    graph::Pathway

    function Metabolism(pathways::Array{Pathway, 1})
        # Build graph by performing all possible Pathway-Pathway merges.
        graph = pathways[1]

        new(pathways, graph)
    end
end


"""
    merge_pathways(pathway1, pathway2)

Create new Pathway from the nodes and edges of two Pathways, but
any shared components are reduced to single copies.  If there are
no shared components, the result will be a disconnected graph containing
full copies of the input Pathways.

# Arguments
- pathway1: Pathway
- pathway2: Pathway

# Returns
- Pathway
"""
function merge_pathways(pathway1, pathway2)
    # Collect nodes and edges into dictionaries for fast lookup.

    # merged_pathway = Pathway(molecules, reactions, reactant_edges, product_edges)

    # return merged_pathway
end
