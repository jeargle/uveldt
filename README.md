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

`Genomes` in &mu;Veldt are simply long strings which can contain `Genes` that code for specific chemical reactions, pores, or transporters.  All `Genes` start with `(` and end at the next occurence of `)`.  For example, the string `AB(CABC)AB)` contains the single `Gene`: `(CABC)`.

`Reactions` specify a set of reactants and a set of products.  They contain element characters as well as the join (`*`) and split (`/`) `Reaction` operators.  `*` means that the elements to either side are separate in the reactants but bonded in the products.  `/` means that they are bonded in the reactants and separate in the products.  To create a `Reaction` string for the inverse of a given `Reaction`, just switch all joins to splits and vice versa.  Multiple joins and/or splits may be defined within a `Reaction` `Gene`.

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

`Genomes` can be built a couple of different ways: through `genome_string()` which creates a new genome string from scratch or by reading a genome from a FASTA file with `read_fasta()`.  Here is an example FASTA file containing two genomes using a chemistry with elements A and B:

    >genome1
    /(*(/B)BBB(/(B(A*A()*)(/B)A)/((B/A)BB*(//B)*A)B(/()B*)A)(*BB(BA*AB(/*/*)*B*/*A/B
    *BA*/(/BA*)AB)(/BB*/AB)A/(BAAA*))AA(/A)*//))A(()(B*)/)B))))A(A/A)(B(////*BA*//A)
    *///B(/A(/()ABAB/*)*ABAB/((B(ABB)/*/A*/)*/(ABA)//*AAABA)*((*())*/)A(B(/AAAA*B)A(
    A*)A*A)AA)*B(//)/A*B)/AB/AB**/*/A*(**)B)((/(A(()B((/BBA)B*A*ABBA((AB*A/**A(/B/A(
    AA)))AA(**AA(**)//)AA*A**(BB///(A(B*B**(B/A)BB)(((A)B((//)A(B)A///B)A*)B*B/***((
    /(**A/*))/AAB))AA(BA)A)B(()B///*ABAB/BA)((A/*AA/(BAB*A(/BB(A)*BA*A/A*/A()A)(//)(
    BA(A)***B*)(B/)/*A)B/)A)B*)/A(*BB)/A(AA()/)B/AB/)*(**(A)(A/**A/)*(BA*/*AB()//B)B
    A/B)B)/B/A)((/*/(/)BB/*AA/B()B)()()*(AA/B))(/(*()(B()A/*(A)(()B/*BB/(**((A*BA*B)
    *(()B/(*))/)AA()B(*//**//(/*/////))/ABB)A/A(A/*A())*(BB)BAB*((B(B//()*)*)/(((ABB
    (*BA))/)*B()***AB))(B()*)B//()(A*A)A)))(A/A)(/BBB*A/B))*(*()*A)B)*(A/(B//A)///A/
    ())**BAA*)BB/)*B*/AB)))(B/(A/(/B/BB))(A*A*)((B/*///)()(ABABAA)///BA*B(((*)*))AA(
    ((*/))*(/*AB/B***(AA(*/)B*A))/*/****/B*/AA/())BB)*A)*ABB/ABA/(A)*AA/B(A*)//()AAB
    A*))*//)/)(//B/*B/B(/AB(A**A(B*((B/*(())
    >genome2
    ((A)A*B//A((/B)*/A(//))/)AB(*/)B)*)(A)/()A(A)/)BB*//(/)A*B)*)B)/A*AAABA)*AAB*((B
    /BB)((*)(*AA*((B*/*(B)(AA()B/A*)*()B)*A)/A(/**BA)BB/)AB)**()*B/A*(//AAA***(/A/B/
    AB*()ABAB/*(///A(BA*BBBBB)*/()/A*(/(*((())B)*B)/)AA/)(**BB/(A(*AA(B)B(A((A(A*BAB
    *A)()AB)(((BB(*))/A)()B*)/)(AABB())()AA/))*AA(()A)(ABB//)B/)(/BB*B)*)A***A)B)/AA
    )B(/*A)*/((A*BABA((B(*B))AA(BAAA*(BB((B/B)AA(A((AB(/)B**/(AA(B()*A)A)))*)A//)B*A
    )*)(A())()AB))A//B))/A*B/(/B/)/B/()(*(*(B(*(/**)((//*(*AB//(B/)**B*)/AB))(B)(ABB
    *AAB(B)B)/AB)/*A/A/B*//)*)A)*/)B)*)(AB)B(A**A*A*//A*BA/BA)*AA((BB(A(/***)/A*)/)A
    )A*B/B*B**B(A*()/*/**A((*B)B/*((A)(*B/B)A(A)B/AA/*(()BB(B/)AA(**/*/B//(*(*)/)(BB
    A((*//A))/B/)/A**ABA*BBBB//))A()AABAA(/(BB/*/BBA(/B)(/B//)(AB*((BBBA*BB)/AA))/((
    *AA*B)(A)A/AB/)))A(AA/)(*(()/(*(/B)*/AAAB)*(/A/A(()/B)*(B*A))B/(B*(BB**((*A)AABA
    )B)//B/(/(A/*(BB/))BA(B*A*/AA/(*(/)BAAB)B(*/B/B)*//*BA)(//(AA(A))*)B//B/B((A/*/B
    B//)*)()(ABA/*A/BA)(B*/B*/*B*)/(/()B(B*B)B*(A*AA//()/*AA(*)*B)B*)B))(*A//(*(*/*)
    (/*((/A/*)*(())/*A//B*)(A*/B)BA))/*(/(A/

You can write one or more genomes to a FASTA file with `write_fasta()`.

When you have a `Genome` loaded, you can use `find_genes()` to discover all of the `Genes` then turn them into their corresponding `Reactions` through calls to `transcribe_gene()` and `parse_reaction()`.


Metabolism
----------

The set of all `Reactions` coded for by a `Genome` is a `Metabolism`.  The structure of a `Metabolism` is a directed graph where nodes are chemicals and reactions, and edges connect reactants to reactions which are also connected to products.  This results in a switching pattern (chemical -> reaction -> chemical -> reaction -> ...) for any path through the metabolic graph.  A chemical can be a product of multiple reactions as well as a reactant for multiple reactions.

This style of graph is a little different from what you might see presented in a biochemistry book where metabolic pathways focus on the "large" molecules and how they change.  These graphs usually have the large molecules as nodes and reactions as edges.  Smaller molecules, like H+, H2O, and CO2 are frequently left out or just shown locally with many copies throughout the graph.  &mu;Veldt handles everything explicitly for completeness and ease of automated downstream analysis and visualization.

`Metabolisms` allow systems biology techniques to be used on &mu;Veldt data.


Evolution
---------

Evolution is biased replication with variation.  In &mu;Veldt, selection is made explicitly using a defined fitness function, and mutations are introduced when the genomes reproduce.  Many types of variation are possible including SNPs, indels, large structural variants, and recombination across multiple genomes.  Support for insertions and deletions means that the genome size is not constant so genetic information can be created or destroyed.  This is uncommon in most genetic algorithms, but it is an important feature of biological evolution.

The genetic algorithm parameters are defined through `SelectionParms` and `MutationParams`.  `SelectionParams` takes a `fitness_function` and one of `select_count`, `select_fraction`, or `fitness_threshold` to determine the final set of passing cells.  Currently, the `fitness_function` must be a function that takes a single `Genome` argument, but this will be extended to cells in the future.  `select_count` limits the passing genomes to a given number, `select_fraction` limits them to a specified fraction of genomes within the population, and `fitness_threshold` only admits those that pass the threshold.

The `MutationParams` consist of mostly intuitive settings for each of the following:

* `snv_rate`: rate of Single Nucleotide Variant; individual character changes.
* `substitution_matrix::SubstitutionMatrix`: frequency matrix for all SNV pairs; used with `snv_rate`.
* `insertion_rate`: rate of multi-character insertions into the genome.
* `insertion_size`: size parameter for the insertion strings.
* `deletion_rate`: rate of multi-character deletions from the genome.
* `deletion_size`: size parameter for the deletion strings.
* `duplication_rate`: rate of multi-character copy/paste to the genome.
* `duplication_size`: size parameter for the duplication strings.
* `inversion_rate`: rate of multi-charcter inversion; string is removed, reversed, and reinserted into the genome.
* `inversion_size`: size parameter for the inversion strings.
* `translocation_rate`: rate of multi-character translocation within the genome; string is removed and then reinserted somewhere else in the genome.
* `translocation_size`: size parameter for the translocation strings.
* `crossing_over::Bool`: whether or not crossing over occurs.

All of the `_rate` and `_size` parameters apply to geometric distributions.

Here's an example evolution parameter YAML file:

    selection_params:
      fitness_function: gene_count
      select_count: 10

    mutation_params:
      snv_rate: 0.01
      insertion_rate: 0.0002
      deletion_rate: 0.0002
      duplication_rate: 0.0001
      inversion_rate: 0.00004
      translocation_rate: 0.00004
      crossing_over: false

And this is a `SubstitutionMatrix` file:

         A    B    C    (    )    *    /
    A    0    5    5    1    1    1    1
    B    5    0    5    1    1    1    1
    C    5    5    0    1    1    1    1
    (    2    2    2    0    1    1    1
    )    2    2    2    1    0    1    1
    *    2    2    2    1    1    0    1
    /    2    2    2    1    1    1    0

Notice how the `Element` characters are included along the top row and down the first column.  Frequencies are read in row:column order, e.g. the frequency of changing B to C is 5 and B to ( is 1.  These frequencies get normalized into probabilities when the file is read in.  You can include the probabilities directly in the matrix, but it can be easier to think about relative frequencies where you start with "1" as the baseline and think of other frequencies as multiples of that.


Lattice
-------

The microbes in &mu;Veldt are embedded in a reaction-diffusion simulation where molecules move randomly throughout the system.  The "world" where these digital microbes live is a 2D or 3D lattice (`Veldt`), and cells live at specific lattice points (`VeldtPoints`).  Across discrete time steps, molecules perform random walks from point to point through the lattice.  If a molecule is located at a point shared by a cell, and the cell is permeable to that molecule, there is some probability that the molecule will enter the cell in the next time step, and vice versa.  Biochemistry encoded by a cell's `Genome` can happen to available reactant molecules within the cell.  Chemistry only occurs within cells.  For now, cells are fixed at their initial lattice locations, and they only interact with the environment though molecular diffusion.

A lot of information is needed to set up a `Veldt`.  Here is example YAML input file that specifies the initial state of a `Veldt`:

    dimensions: [3, 4]

    chemistry: ./chemistries/chemistry1.yml

    molecules:
      - name: AAA
        distribution:
        - type: even
        - count: 10
      - name: BBB
        distribution:
        - type: even
        - count: 30

    genomes:
      - ./genomes/genomes1.fasta

    cells:
      - genome: genome1
        location: [1, 1]
        molecules:
        - name: AAA
          count: 11
        - name: BBB
          count: 22
      - genome: genome2
        location: [2, 2]
        molecules:
        - name: AA
          count: 33
        - name: AB
          count: 44

In this `Veldt`, the lattice structure is defined by `dimensions`, which specifies a 2D lattice with 3 rows and 4 columns (12 `VeldtPoints` total).  To make a 3D lattice, you just need one more number for the "depth" dimension.  The `chemistry` line points to a YAML input file for the `Chemistry`.  The `molecules` section sets AAA and BBB as evenly distributed across the lattice with AAA having 10 molecules at each point and BBB 30.  The `genomes` section points to a FASTA file containing the "genome1" and "genome2" `Genomes`.  More than one FASTA file could be included here.  The `cells` section defines two cells: the first with `Genome` genome1 at `VeldtPoint` (1, 1) containing 11 AAA molecules and 22 BBB molecules, and the secont with `Genome` genome2 at `VeldtPoint` (2, 2) containing 33 AA molecules and 44 AB molecules.  Keeping `Chemistries` and `Genomes` separate in their own files makes it easy to reuse them and limit the focus of the `Veldt` file to the lattice and its contents.  The molecules and `Genomes` are constrained to elements present in the provided `Chemistry` so this `Chemistry` cleary includes elements A and B.


How to Run
----------

&mu;Veldt has a commandline script `uveldtsim.jl` that reads a setup file and runs the specified simulation.  For example, to run the 20-step test simulation in `simulation1.yml`, navigate to `uveldt/bin` and run

```
> ./uveldtsim.jl ../test/simulations/simulation1.yml
```

Data for each step of the simulation will be printed to stdout.


Dependencies
------------

* ArgParse
* DataStructures
* Distributions
* LinearAlgebra
* Printf
* Random
* UUIDs
* YAML
