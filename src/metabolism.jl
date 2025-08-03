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

        Metabolism(nodes, edges)
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
- `Metabolism`: merged Metabolism
"""
function merge_metabolisms(metabolism1, metabolism2)
    nodes = union(metabolism1.nodes, metabolism2.nodes)
    edges = union(metabolism1.edges, metabolism2.edges)

    return Metabolism(nodes, edges)
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
function get_connected_submetabolism(metabolism::Metabolism, node::Union{Molecule, Reaction})
    final_nodes = Set{Union{Molecule, Reaction}}()
    final_edges = Set{Tuple{Union{Molecule, Reaction}, Union{Molecule, Reaction}}}()

    next_nodes = Set{Union{Molecule, Reaction}}([node])

    while length(next_nodes) > 0
        temp_nodes = next_nodes
        next_nodes = Set{Union{Molecule, Reaction}}()

        for temp_node in temp_nodes
            if !(temp_node in final_nodes)
                children = metabolism.children[temp_node]
                parents = metabolism.parents[temp_node]
                union!(next_nodes, children)
                union!(next_nodes, parents)

                push!(final_nodes, temp_node)

                for child in children
                    push!(final_edges, (temp_node, child))
                end

                for parent in parents
                    push!(final_edges, (parent, temp_node))
                end
            end
        end
    end

    return Metabolism(final_nodes, final_edges)
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
    submetabolisms = Array{Metabolism, 1}()
    used_nodes = Set{Union{Molecule, Reaction}}()

    for node in metabolism.nodes
        if !(node in used_nodes)
            submet = get_connected_submetabolism(metabolism, node)
            union!(used_nodes, submet.nodes)
            push!(submetabolisms, submet)
        end
    end

    return submetabolisms
end
