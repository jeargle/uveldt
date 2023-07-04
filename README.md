&mu;Veldt
=========

Evolution and biochemistry of digital microbes


General
-------

The basic idea for this project is to simulate microbial evolution using simplified biochemistry and genetics.

All chemistry is based on string splitting and concatenation.  Genomes code directly for operators (analogous to proteins) that accept reactant strings and generate product strings.  Chemical energy is used and produced by these reactions so that there is some measure of how fit an organism is.  Each organism is basically a cell that contains a genome and a set of molecules.  Molecules diffuse on a grid, enter and leave cells, and are acted upon by the cytosolic proteins.

Organisms are judged based on their energy pools at the end of a generation, then the winners reproduce with mutation, and the next generation begins.  Mutation can include structural variants, like large scale insertions or deletions, as well as SNPs.  Genomes are long strings that can code for multiple proteins, but they may also contain pseudogenes and runs of noncoding DNA.  This allows genomes to gain and lose functionality over generations and makes them more robust to mutation.


Chemistry
---------

The chemistry of &mu;Veldt is based around simple string operations.  Elements are represented as single characters, molecules are 1-dimensional strings of elements, and bonds connect pairs of elements together.  All chemical reactions involve making and/or breaking bonds between the elements in molecules.  For example, atoms of elements `A` and `B` can be connected to make the molecule `AB`.

A `Chemistry` is a set of `Element` definitions and parameters for pairwise bonds between them.  These definitions can be set up in YAML files with a form like:

    elements:
      - name: A
        mass: 1
      - name: B
        mass: 2
      - name: C
        mass: 3

    bonds:
      - elements:
        - A
        - B
        energy_change: -5.5
        transition_energy: 1.5
      - elements:
        - A
        - C
        energy_change: -6.5
        transition_energy: 2.5
      - elements:
        - B
        - C
        energy_change: -7.5
        transition_energy: 3.5

Each `Element` has a mass, which defaults to 1.  Bonds have both an `energy_change`, which is the change in energy after bond formation, and a `transition_energy` needed up front to start the bond creation process.  These energies are tracked within each cell during simulations.  A complete `Chemistry` requires bond definitions for every pair of `Elements`.


Genome
------

Genomes in &mu;Veldt are simply long strings which can contain genes that code for specific chemical reactions, pores, or transporters.  All genes start with `(` and end at the next occurence of `)`.  For example, the string `AB(CABC)AB)` contains the single gene: `(CABC)`.

Reactions specify a set of reactants and a set of products.  They contain element characters as well as the join (`*`) and split (`/`) reaction operators.  `*` means that the elements to either side are separate in the reactants but bonded in the products.  `/` means that they are bonded in the reactants and separate in the products.  To create a reaction string for the inverse of a given reaction, just switch all joins to splits and vice versa.  Multiple joins and/or splits may be defined within a reaction gene.

Examples:

* (A\*B): A + B -> AB
* (A/B): AB -> A + B
* (AB\*C): AB + C -> ABC
* (AB/C): ABC -> AB + C
* (A\*B/C): A + BC -> AB + C
* (A\*B\*C): A + B + C -> ABC
* (A\*A\*A): A + A + A -> AAA

Pores are specified by single molecule strings with no join or split operators.  If a pore is active, molecules of that type are allowed to diffuse freely into and out of the cell.

Examples:

* (A): A pore
* (AB): AB pore
* (AAA): AAA pore

Transporters are similar to pores except that they use energy to move molecules across the cell membrane.  A transporter string must start with a join or split operator followed by a single molecule string.

Examples:

* (\*A): A(out) -> A(in)
* (\*AB): AB(out) -> AB(in)
* (/A): A(in) -> A(out)
* (/AAA): AAA(in) -> AAA(out)


Evolution
---------

Evolution is biased replication with variation.  In &mu;Veldt, selection is made explicitly using a defined fitness function, and mutations are introduced when the genomes reproduce.  Many types of variation are possible including SNPs, indels, large structural variants, and recombination across multiple genomes.  Support for insertions and deletions means that the genome size is not constant so genetic information can be created or destroyed.  This is uncommon in most genetic algorithms, but it is an important feature of biological evolution.

The genetic algorithm parameters are defined through `SelectionParms` and `MutationParams`.  `SelectionParams` takes a `fitness_function` and one of `select_count`, `select_fraction`, or `fitness_threshold` to determine the final set of passing cells.  Currently, the `fitness_function` must be a function that takes a single `Genome` argument, but this will be extended to cells in the future.  `select_count` limits the passing genomes to a given number, `select_fraction` limits them to a specified fraction of genomes within the population, and `fitness_threshold` only admits those that pass the threshold.

The `MutationParams` consist of mostly intuitive settings for each of the following:

* `snv_rate`
* `substitution_matrix::SubstitutionMatrix`
* `insertion_rate`
* `insertion_size`
* `deletion_rate`
* `deletion_size`
* `duplication_rate`
* `duplication_size`
* `inversion_rate`
* `inversion_size`
* `translocation_rate`
* `translocation_size`
* `crossing_over::Bool`: whether or not crossing over occurs

All of the `_rate` and `_size` parameters apply to geometric distributions.


Lattice
-------

The microbes in &mu;Veldt are embedded in a reaction-diffusion simulation where molecules move randomly throughout the system.  The "world" where these digital microbes live is a 2D or 3D lattice (`Veldt`), and cells live at specific lattice points (`VeldtPoints`).  Across discrete time steps, molecules perform random walks from point to point through the lattice.  If a molecule is located at a point shared by a cell, and the cell is permeable to that molecule, there is some probability that the molecule will enter the cell in the next time step, and vice versa.  Biochemistry encoded by a cell's genome can happen to available reactant molecules within the cell.  Chemistry only occurs within cells.  For now, cells are fixed at their initial lattice locations, and they only interact with the environment though molecular diffusion.


Dependencies
------------

* Distributions
* Printf
* Random
* YAML
