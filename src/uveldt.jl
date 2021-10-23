# John Eargle (mailto: jeargle at gmail.com)
# 2018-2021
# uveldt

module uveldt

using Distributions
using Printf
using Random
using YAML

using LinearAlgebra
using UUIDs

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

# Evolution
include("evolution.jl")
export MutationParams, SubstitutionMatrix
export read_substitution_matrix
export add_snvs, add_insertions, remove_deletions, add_duplications
export add_inversions, add_translocations, cross_over, mutate

# Simulation
include("simulation.jl")
export Cell, VeldtPoint, Veldt
export init_molecules, add_cell, get_neighbors, setup_veldt, simulate

end
