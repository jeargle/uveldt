# John Eargle (mailto: jeargle at gmail.com)
# uveldt

module uveldt

using DataStructures
using Distributions
using Printf
using Random
using YAML

using LinearAlgebra
using UUIDs

# utils
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end

# Chemistry
include("chemistry.jl")
export Element, ElementTable, Bond, BondTable, Chemistry, Molecule
export ReactionType, Reaction

export get_bond, mass, setup_chemistry

# Genome
include("genome.jl")
export Gene, Genome
export parse_reaction, transcribe_gene, genome_string, find_genes
export is_pseudogene, read_fasta, write_fasta

# Metabolism
include("metabolism.jl")
export Metabolism
export merge_pathways

# Evolution
include("evolution.jl")
export PhyloNode, PhyloEdge, Phylogeny, SelectionParams, MutationParams, SubstitutionMatrix
export select_genomes, select_cells
export gene_count, genome_length, fitness_functions
export read_substitution_matrix, read_evolution_params
export add_snvs, add_insertions, remove_deletions, add_duplications
export add_inversions, add_translocations, cross_over, mutate

# Simulation
include("simulation.jl")
export Cell, VeldtPoint, Veldt
export init_molecules, add_cell, get_neighbors, setup_veldt
export simulate, setup_simulation

end
