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

Genomes in &mu;Veldt are simply long strings which can contain genes that code for specific chemical reactions.  These reactions specify a set of reactants and a set of products.  Genes start with "(" and end with ")".  They contain element characters as well as the reaction operators "\*" and "/".  "\*" means that the elements to either side are separate in the reactants but bonded in the products.  "/" means that they are bonded in the reactants and separate in the products.

Examples:

* (A\*B): A + B -> AB
* (A/B): AB -> A + B
* (A\*B/C): A + BC -> AB + C
* (A\*B\*C): A + B + C -> ABC
* (A\*A\*A): A + A + A -> AAA


Evolution
---------

Evolution is biased replication with variation.  Selection is made using some fitness function and mutations are introduces when the genomes reproduce.  Many types of variation are possible including SNPs, indels, large structural variants, and recombination across multiple genomes.
