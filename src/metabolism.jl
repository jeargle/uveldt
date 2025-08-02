# John Eargle (mailto: jeargle at gmail.com)
# uveldt.metabolism


"""
Graph with nodes representing Molecules and Reactions.
"""
struct Metabolism
    nodes::Set{Union{Molecule, Reaction}}
    edges::Set{Tuple{Union{Molecule, Reaction}, Union{Molecule, Reaction}}}
    children::Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}
    parents::Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}

    function Metabolism(nodes, edges, children, parents)
        # Used for taking subsets of existing Metabolisms.
        new(nodes, edges, children, parents)
    end

    function Metabolism(nodes::Set{Union{Molecule, Reaction}}, edges::Set{Tuple{Union{Molecule, Reaction}, Union{Molecule, Reaction}}})
        children = Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}()
        parents = Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}()

        for node in nodes
            children[node] = Set{Union{Molecule, Reaction}}()
            parents[node] = Set{Union{Molecule, Reaction}}()
        end

        for (parent, child) in edges
            push!(children[parent], child)
            push!(parents[child], parent)
        end

        new(nodes, edges, children, parents)
    end

    function Metabolism(reactions::Array{Reaction, 1})
        nodes = Set{Union{Molecule, Reaction}}()
        edges = Set{Tuple{Union{Molecule, Reaction}, Union{Molecule, Reaction}}}()

        for reaction in reactions
            push!(nodes, reaction)
            for reactant in reaction.reactants
                push!(nodes, reactant)
                push!(edges, (reactant, reaction))
            end

            for product in reaction.products
                push!(nodes, product)
                push!(edges, (reaction, product))
            end
        end

        new(nodes, edges)
    end
end


"""
    get_source_nodes(metabolism)

Retrieve a set of nodes that have no parents.

# Arguments
- `metabolism`: Metabolism

# Returns
- `Set{Union{Molecule, Reaction}}`:
"""
function get_source_nodes(metabolism)
    source_nodes = Set{Union{Molecule, Reaction}}()

    for node in metabolism.nodes
        if !haskey(metabolism.parents, node)
            push!(source_nodes, node)
        end
    end

    return source_nodes
end


"""
    get_terminal_nodes(metabolism)

Retrieve a set of nodes that have no children.

# Arguments
- `metabolism`: Metabolism

# Returns
- `Set{Union{Molecule, Reaction}}`:
"""
function get_terminal_nodes(metabolism)
    terminal_nodes = Set{Union{Molecule, Reaction}}()

    for node in metabolism.nodes
        if !haskey(metabolism.children, node)
            push!(terminal_nodes, node)
        end
    end

    return terminal_nodes
end


"""
    merge_metabolisms(metabolism1, metabolism2)

Create new Metabolism from the nodes and edges of two Metabolisms, but
any shared components are reduced to single copies.  If there are
no shared components, the result will be a disconnected graph containing
full copies of the input Metabolisms.

# Arguments
- `metabolism1`: Metabolism
- `metabolism2`: Metabolism

# Returns
- `Metabolism`:
"""
function merge_metabolisms(metabolism1, metabolism2)
    # Collect nodes and edges into dictionaries for fast lookup.

    # merged_metabolism = Metabolism(molecules, reactions, reactant_edges, product_edges)

    # return merged_metabolism
end


"""
    get_connected_metabolism(metabolism)

Create a Metabolism from the connected subgraph of an input Metabolism and one
of its nodes.

# Arguments
- `metabolism`: Metabolism
- `node`: Union{Molecule, Reaction}

# Returns
- `Metabolism`:
"""
function get_connected_submetabolism(metabolism, node)
    # Collect nodes and edges into dictionaries for fast lookup.

    next_nodes = Set{Union{Molecule, Reaction}}()
    union!(next_nodes, metabolism.children[node])
    union!(next_nodes, metabolism.parents[node])
    # submetabolism = Metabolism(molecules, reactions, reactant_edges, product_edges)

    # return submetabolism
end


"""
    get_connected_metabolisms(metabolism)

Create a list of Metabolisms from the connected subgraphs of an input
Metabolism.

# Arguments
- `metabolism`: Metabolism

# Returns
- `Array{Metabolism, 1}`:
"""
function get_connected_submetabolisms(metabolism)
    # Collect nodes and edges into dictionaries for fast lookup.

    submetabolisms = Array{Metabolism, 1}()
    used_nodes = Set{Union{Molecule, Reaction}}()

    for node in metabolism.nodes
        if node not in used_nodes
            submet = get_connected_submetebolism(metabolism, node)
            union!(used_nodes, submet.nodes)
            push!(submetabolisms, submet)
        end
    end

    return submetabolisms
end
