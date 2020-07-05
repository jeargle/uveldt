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

The chemistry of &mu;Veldt is based around simple string operations.  Elements are represented as single characters, molecules are 1-dimensional strings of elements, and bonds connect pairs of elements together.  All chemical reactions involve making and/or breaking bonds between the elements in molecules.

For example, atoms of elements "A" and "B" can be connected to make the molecule "AB".


Genome
------

Genomes in &mu;Veldt are simply long strings which can contain genes that code for specific chemical reactions, pores, or transporters.

Reactions specify a set of reactants and a set of products.  Genes start with "(" and end with ")".  They contain element characters as well as the reaction operators "\*" and "/".  "\*" means that the elements to either side are separate in the reactants but bonded in the products.  "/" means that they are bonded in the reactants and separate in the products.  To create a reaction string for the inverse of a given reaction, just switch all joins to splits and vice versa.

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


Lattice
-------

The microbes in &mu;Veldt are embedded in a reaction-diffusion simulation where molecules move randomly throughout the system.  The "world" where these digital microbes live is a 2D or 3D lattice (Veldt), and cells live at specific lattice points (VeldtPoints).  Across discrete time steps, molecules perform random walks from point to point through the lattice.  If a molecule is locate at a point shared by a cell, and the cell is permeable to that molecule, there is some probability that the molecule will enter the cell in the next time step.  Biochemistry encoded by a cell's genome can happen to available reactant molecules within the cell.  Chemistry only occurs within cells.  For now, cells are fixed at their initial lattice locations so they only interact with the environment though molecular diffusion.


Dependencies
------------

* Distributions
* Printf
* Random
* YAML
