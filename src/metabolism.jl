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

    function Metabolism(reactions::Array{Reaction, 1})
        nodes = Set{Union{Molecule, Reaction}}()
        edges = Set{Tuple{Union{Molecule, Reaction}, Union{Molecule, Reaction}}}()
        children = Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}()
        parents = Dict{Union{Molecule, Reaction}, Set{Union{Molecule, Reaction}}}()

        for reaction in reactions
            push!(nodes, reaction)
            for reactant in reaction.reactants
                push!(nodes, reactant)
                push!(edges, (reactant, reaction))
                if haskey(children, reactant)
                    push!(children[reactant], reaction)
                else
                    children[reactant] = Set([reaction])
                end
                if haskey(parents, reaction)
                    push!(parents[reaction], reactant)
                else
                    parents[reaction] = Set([reactant])
                end
            end
            for product in reaction.products
                push!(nodes, product)
                push!(edges, (reaction, product))
                if haskey(children, reaction)
                    push!(children[reaction], product)
                else
                    children[reaction] = Set([product])
                end
                if haskey(parents, product)
                    push!(parents[product], reaction)
                else
                    parents[product] = Set([reaction])
                end
            end
        end
        new(nodes, edges, children, parents)
    end
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
