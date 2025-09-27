# John Eargle (mailto: jeargle at gmail.com)
# test
#
# To build sysimage boom.so from uveldt/test:
#   using PackageCompiler
#   create_sysimage([:Distributions, :YAML, :Printf, :Random], sysimage_path="../boom.so", precompile_execution_file="so_builder.jl")
#
# To run from uveldt/test:
#   julia --project=.. -J../boom.so test.jl


using DataStructures
using Printf
using Test
using uveldt


function print_test_header(test_name)
    border = repeat("*", length(test_name) + 4)
    println(border)
    println("* ", test_name, " *")
    println(border)
end

function create_elements(el_count)
    println("  * Create Elements")

    alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    elements = []

    for i in 1:el_count
        push!(elements, Element(alphabet[i], i))
    end

    return elements
end

function create_bonds(element_table, energies)
    println("  * Create Bonds")

    bonds = []

    energy_idx = 1
    elements = sort([i for i in keys(element_table.elements)])
    el_count = length(elements)
    for i in 1:el_count
        el1 = element_table.elements[elements[i]]
        for j in i:el_count
            el2 = element_table.elements[elements[j]]
            energy1, energy2 = energies[energy_idx]
            push!(bonds, Bond(el1, el2, energy1, energy2))
            energy_idx += 1
        end
    end

    return bonds
end

function test_element()
    print_test_header("Element")

    el1, el2, el3 = create_elements(3)

    println("el1: ", el1)
    println("el2: ", el2)
    println("el3: ", el3)
    println()

    @test el1.name == 'A'
    @test el2.name == 'B'
    @test el3.name == 'C'

    println("mass(el1): ", mass(el1))
    println("mass(el2): ", mass(el2))
    println("mass(el3): ", mass(el3))
    println()

    @test mass(el1) == 1
    @test mass(el2) == 2
    @test mass(el3) == 3
end

function test_element_table()
    print_test_header("ElementTable")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    println("el_table1: ", el_table1)
    println("el_table1.elements['A']: ", el_table1.elements['A'])
    println("el_table1.elements['B']: ", el_table1.elements['B'])
    println()

    @test el_table1.elements['A'] == elements[1]
    @test el_table1.elements['B'] == elements[2]
end

function test_bond()
    print_test_header("Bond")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bond1, bond2, bond3 = create_bonds(el_table1, energies1)
    println("bond1: ", bond1)
    println("bond2: ", bond2)
    println("bond3: ", bond3)
    println()

    @test bond1.energy_change == -5.5
    @test bond1.transition_energy == 1.5
    @test bond2.energy_change == -4.5
    @test bond2.transition_energy == 0.5
    @test bond3.energy_change == 3.6
    @test bond3.transition_energy == 3.4
end

function test_bond_table()
    print_test_header("BondTable")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)
    println("b_table1: ", b_table1)

    println("  * get Bonds")
    println("get_bond(b_table1, 'A', 'A'): ", get_bond(b_table1, 'A', 'A'))
    println("get_bond(b_table1, 'B', 'B'): ", get_bond(b_table1, 'B', 'B'))
    println("get_bond(b_table1, 'A', 'B'): ", get_bond(b_table1, 'A', 'B'))
    println("get_bond(b_table1, 'B', 'A'): ", get_bond(b_table1, 'B', 'A'))
    println()

    @test get_bond(b_table1, 'A', 'A') == bonds[1]
    @test get_bond(b_table1, 'A', 'B') == bonds[2]
    @test get_bond(b_table1, 'B', 'A') == bonds[2]
    @test get_bond(b_table1, 'B', 'B') == bonds[3]
end

function test_chemistry()
    print_test_header("Chemistry")

    elements1 = create_elements(2)
    el_table1 = ElementTable(elements1)

    elements2 = create_elements(3)
    el_table2 = ElementTable(elements2)

    println("  * add Bonds")
    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (-4.5, 0.5), (3.6, 3.4), (3.6, 3.4), (3.6, 3.4)]
    bonds = create_bonds(el_table2, energies1)
    bond1, bond2, bond3, bond4, bond5, bond6 = bonds

    b_table1 = BondTable([bond1, bond2, bond4])
    b_table2 = BondTable(bonds)

    println("  * create Chemistries")
    chem1 = Chemistry(el_table1, b_table1)
    chem2 = Chemistry(el_table2, b_table2)
    println("chem1: ", chem1)
    println("chem2: ", chem2)
    println()

    @test chem1.element_table == el_table1
    @test chem1.bond_table == b_table1
    @test chem2.element_table == el_table2
    @test chem2.bond_table == b_table2
end

function test_chemistry_setup()
    print_test_header("Chemistry Setup")

    chem1 = setup_chemistry("./chemistries/chemistry1.yml")
    chem2 = setup_chemistry("./chemistries/chemistry2.yml")

    println("chem1: ", chem1)
    println("chem2: ", chem2)

    element1 = chem1.element_table.elements['A']
    element2 = chem1.element_table.elements['B']
    element3 = chem1.element_table.elements['C']
    @test element1.name == 'A'
    @test element1.mass == 1
    @test element2.name == 'B'
    @test element2.mass == 2
    @test element3.name == 'C'
    @test element3.mass == 3

    bond1 = get_bond(chem1.bond_table, 'A', 'A')
    bond2 = get_bond(chem1.bond_table, 'A', 'B')
    bond3 = get_bond(chem1.bond_table, 'A', 'C')
    @test bond1.energy_change == -5.5
    @test bond1.transition_energy == 1.5
    @test bond2.energy_change == -5.5
    @test bond2.transition_energy == 1.5
    @test bond3.energy_change == -6.5
    @test bond3.transition_energy == 2.5

    element1 = chem2.element_table.elements['H']
    element2 = chem2.element_table.elements['C']
    element3 = chem2.element_table.elements['O']
    @test element1.name == 'H'
    @test element1.mass == 1
    @test element2.name == 'C'
    @test element2.mass == 12
    @test element3.name == 'O'
    @test element3.mass == 16

    bond1 = get_bond(chem2.bond_table, 'H', 'C')
    bond2 = get_bond(chem2.bond_table, 'H', 'O')
    bond3 = get_bond(chem2.bond_table, 'C', 'O')
    @test bond1.energy_change == -4.5
    @test bond1.transition_energy == 1.1
    @test bond2.energy_change == -6.5
    @test bond2.transition_energy == 1.3
    @test bond3.energy_change == -8.5
    @test bond3.transition_energy == 1.5

end

function test_molecule()
    print_test_header("Molecule")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    mol1 = Molecule("mol1", "AAA", chem1)
    mol2 = Molecule("mol2", "BBB", chem1)
    mol3 = Molecule("mol3", "ABA", chem1)
    mol4 = Molecule("BAB", chem1)
    mol5 = Molecule("AAB", chem1)
    mol6 = Molecule("BBA", chem1)
    println("mol1: ", mol1)
    println("mol2: ", mol2)
    println("mol3: ", mol3)
    println("mol4: ", mol4)
    println("mol5: ", mol5)
    println("mol6: ", mol6)
    println()

    @test mol1.name == "mol1"
    @test mol2.name == "mol2"
    @test mol3.name == "mol3"
    @test mol4.name == "BAB"
    @test mol5.name == "AAB"
    @test mol6.name == "BBA"

    println("mass(mol1): ", mass(mol1))
    println("mass(mol2): ", mass(mol2))
    println("mass(mol3): ", mass(mol3))
    println("mass(mol4): ", mass(mol4))
    println("mass(mol5): ", mass(mol5))
    println("mass(mol6): ", mass(mol6))
    println()

    @test mass(mol1) == 3
    @test mass(mol2) == 6
    @test mass(mol3) == 4
    @test mass(mol4) == 5
    @test mass(mol5) == 4
    @test mass(mol6) == 5

end

function test_reaction()
    print_test_header("Reaction")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    println("  * Reactions")
    reactions = [Reaction(["A", "A"], ["AA"], chem1),
                 Reaction(["B", "B"], ["BB"], chem1),
                 Reaction(["AB", "BA"], ["ABBA"], chem1),
                 Reaction(["AA"], ["A", "A"], chem1),
                 Reaction(["BB"], ["B", "B"], chem1),
                 Reaction(["ABBA"], ["AB", "BA"], chem1)]
    for (i, reaction) in enumerate(reactions)
        @printf "reaction%d: %s\n" i reaction
    end

    println("  * Pores")
    reactions = [Reaction("A", chem1),
                 Reaction("BB", chem1),
                 Reaction("ABBA", chem1)]
    for (i, reaction) in enumerate(reactions)
        @printf "pore%d: %s\n" i reaction
    end

    println("  * Transporters")
    reactions = [Reaction("A", chem1, transport_in, -3.0),
                 Reaction("BB", chem1, transport_in, -4.0),
                 Reaction("ABBA", chem1, transport_in, -5.0),
                 Reaction("A", chem1, transport_out, -3.0),
                 Reaction("BB", chem1, transport_out, -4.0),
                 Reaction("ABBA", chem1, transport_out, -5.0)]
    for (i, reaction) in enumerate(reactions)
        @printf "pore%d: %s\n" i reaction
    end

    println()

end

function test_gene()
    print_test_header("Gene")

    elements = create_elements(3)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (-4.5, 0.5), (3.6, 3.4), (3.6, 3.4), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    genes = [Gene(1, "(A*B)", chem1),
             Gene(21, "(A/B)", chem1),
             Gene(45, "(A*B/C)", chem1),
             Gene(132, "(A*B*C)", chem1),
             Gene(154, "(A*//*A)", chem1),
             Gene(174, "(A/*A*/A)", chem1),
             Gene(194, "(A)", chem1),
             Gene(214, "(AB)", chem1),
             Gene(234, "(*AB)", chem1),
             Gene(254, "(/BA)", chem1),
             Gene(274, "(*/AAA/*/*/)", chem1),
             Gene(294, "(/*B*B/B*/*/*)", chem1),
             Gene(314, "()", chem1)]

    for i in 1:length(genes)
        @printf "gene%d: %s\n" i genes[i]
    end

    println("  * Gene transcripts")
    reactions = []
    for i in 1:length(genes)
        @printf "gene%d.transcript: %s\n" i genes[i].transcript
        reaction_type, reactants, products = parse_reaction(genes[i].transcript)
        println("reaction type: ", reaction_type)
        println("reactants: ", reactants)
        println("products: ", products)
        if is_pseudogene(genes[i])
            println("  *** pseudogene: " * genes[i].string)
        elseif reaction_type == reaction
            push!(reactions, Reaction(reactants, products, chem1))
        elseif reaction_type == pore
            push!(reactions, Reaction(reactants[1], chem1))
        else
            push!(reactions, Reaction(reactants[1], chem1, reaction_type, -3.0))
        end
    end

    println("  * Reactions")
    for i in 1:length(reactions)
        @printf "reaction%d: %s\n" i reactions[i]
    end

    println()
end

function print_genes(genes)
    for gene in genes
        if gene.strand == positive
            strand_str = "+"
        else
            strand_str = "-"
        end

        @printf "%s%d %s\n" strand_str gene.location gene.string
    end
end

function test_genome()
    print_test_header("Genome")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    genome_str1 = genome_string(1000, chem1)
    println("genome_str1: ", genome_str1)

    genome1 = Genome("genome1", genome_str1, chem1)
    println("genome1: ", genome1)
    write_fasta(genome1, "genome1.fasta")

    genome_info = read_fasta("genome1.fasta")
    println("length(genome_info): ", length(genome_info))
    genome_str2 = genome_info[1][3]
    println("genome_str2: ", genome_str2)
    println("name1 == name2: ", genome_info[1][1] == "genome1")
    println("genome_str1 == genome_str2: ", genome_str1 == genome_str2)

    genes = find_genes(genome1)
    genes2 = find_genes(genome1, include_reverse=true)

    println("  * find Genes")
    print_genes(genes)

    println("  * find Genes 2")
    print_genes(genes2)

    println("  * Gene transcripts")
    pseudogene_count = 0
    reaction_count = 0
    pore_count = 0
    transport_in_count = 0
    transport_out_count = 0
    for i in 1:length(genes)
        @printf "gene%d.transcript: %s\n" i genes[i].transcript
        if is_pseudogene(genes[i])
            println("  *** pseudogene: " * genes[i].string)
            pseudogene_count += 1
        else
            reaction_type, reactants, products = parse_reaction(genes[i].transcript)
            println("  reaction type: ", reaction_type)
            println("  reactants: ", reactants)
            println("  products: ", products)
            if reaction_type == reaction
                reaction_count += 1
            elseif reaction_type == pore
                pore_count += 1
                @test length(reactants) == 1
            elseif reaction_type == transport_in
                transport_in_count += 1
                @test length(reactants) == 1
            elseif reaction_type == transport_out
                transport_out_count += 1
                @test length(reactants) == 1
            end
        end
    end

    @printf "%d good genes\n" length(genes)-pseudogene_count
    @printf "  %d reactions\n" reaction_count
    @printf "  %d pores\n" pore_count
    @printf "  %d transport in\n" transport_in_count
    @printf "  %d transport out\n" transport_out_count
    @printf "%d pseudogenes\n" pseudogene_count
    @printf "%d total\n" length(genes)

    genome_str2 = genome_string(500, chem1)
    println("genome_str2: ", genome_str2)
    genome2 = Genome("genome2", genome_str2, chem1)
    println("genome2: ", genome2)

    println()
    println("  * Add SNVs")
    genome3 = add_snvs(genome1, 0.01)
    println("genome3: ", genome3)
    println()
    println("  * Add insertions")
    genome4 = add_insertions(genome1, 0.002)
    println("genome4: ", genome4)
    println()
    println("  * Remove deletions")
    genome5 = remove_deletions(genome1, 0.002)
    println("genome5: ", genome5)
    println()
    println("  * Add duplications")
    genome6 = add_duplications(genome1, 0.002, size_param=0.02)
    println("genome6: ", genome6)
    println()
    println("  * Add inversions")
    genome7 = add_inversions(genome1, 0.002, size_param=0.02)
    println("genome7: ", genome7)
    println()
    println("  * Add translocations")
    genome8 = add_translocations(genome1, 0.002, size_param=0.02)
    println("genome8: ", genome8)

    println()
    println("  * Cross over")
    genome9, genome10 = cross_over(genome1, genome2)
    println("genome9: ", genome9)
    println("genome10: ", genome10)

    println()
end

function test_cell()
    print_test_header("Cell")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)
    mol1 = Molecule("AAA", chem1)
    mol2 = Molecule("BBB", chem1)

    genome_str = genome_string(200, chem1)
    println("genome_string: ", genome_str)

    genome1 = Genome("genome1", genome_str, chem1)
    println("genome1: ", genome1)

    println("  * Cells")
    cell1 = Cell(genome1)
    println("cell1: ", cell1)

    println("  * add Molecule counts")
    cell1.molecule_counts[1][mol1] = 33
    cell1.molecule_counts[1][mol2] = 44
    println("cell1: ", cell1)

    println()
end

function test_veldt_point()
    print_test_header("VeldtPoint")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)
    mol1 = Molecule("AAA", chem1)
    mol2 = Molecule("BBB", chem1)

    genome_str = genome_string(200, chem1)
    genome1 = Genome("genome1", genome_str, chem1)
    cell1 = Cell(genome1)
    cell1.molecule_counts[1][mol1] = 33
    cell1.molecule_counts[1][mol2] = 44

    println("  * VeldtPoints")
    veldt_pt1 = VeldtPoint()
    veldt_pt2 = VeldtPoint()
    veldt_pt3 = VeldtPoint(cell=cell1)
    println("veldt_pt1: ", veldt_pt1)
    println("veldt_pt2: ", veldt_pt2)
    println("veldt_pt3: ", veldt_pt3)

    println("  * add Molecule counts")
    veldt_pt1.molecule_counts[1][mol1] = 300
    veldt_pt1.molecule_counts[1][mol2] = 400
    veldt_pt2.molecule_counts[1][mol1] = 330
    veldt_pt2.molecule_counts[1][mol2] = 440
    veldt_pt2.molecule_counts[1][mol1] = 333
    veldt_pt2.molecule_counts[1][mol2] = 444
    println("veldt_pt1: ", veldt_pt1)
    println("veldt_pt2: ", veldt_pt2)
    println("veldt_pt3: ", veldt_pt3)

    println()
end

function test_veldt()
    print_test_header("Veldt")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.6, 3.4)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)
    mol1 = Molecule("AAA", chem1)
    mol2 = Molecule("BBB", chem1)
    mol3 = Molecule("AA", chem1)
    mol4 = Molecule("BB", chem1)

    genome_str1 = genome_string(200, chem1)
    genome1 = Genome("genome1", genome_str1, chem1)
    genome_str2 = genome_string(200, chem1)
    genome2 = Genome("genome2", genome_str2, chem1)

    cell1 = Cell(genome1)
    cell1.molecule_counts[1][mol1] = 11
    cell1.molecule_counts[1][mol2] = 22
    cell2 = Cell(genome2)
    cell2.molecule_counts[1][mol3] = 33
    cell2.molecule_counts[1][mol4] = 44

    println("  * Veldt")
    veldt1 = Veldt([3, 4])
    println("veldt1: ", veldt1)

    println("  * initialize Molecule counts")

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 10
    dd1[mol2] = 10
    init_molecules(veldt1, [1, 1], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 10
    dd1[mol2] = 20
    init_molecules(veldt1, [1, 2], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 10
    dd1[mol2] = 30
    init_molecules(veldt1, [1, 3], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 20
    dd1[mol2] = 10
    init_molecules(veldt1, [2, 1], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 20
    dd1[mol2] = 20
    init_molecules(veldt1, [2, 2], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 20
    dd1[mol2] = 30
    init_molecules(veldt1, [2, 3], dd1)

    dd1 = DefaultDict{Molecule, Int64}(0)
    dd1[mol1] = 30
    dd1[mol2] = 30
    init_molecules(veldt1, [3, 3], dd1)

    println("veldt1: ", veldt1)

    println("  * add Cells")
    add_cell(veldt1, [1, 1], cell1)
    add_cell(veldt1, [2, 2], cell2)
    println("veldt1: ", veldt1)

    println("  * get neighboring VeldtPoints")
    vps = get_neighbors(veldt1, [1, 1])
    println("  [1, 1]")
    println("-x: ", vps[1])
    println("+x: ", vps[2])
    println("-y: ", vps[3])
    println("+y: ", vps[4])

    vps = get_neighbors(veldt1, [2, 2])
    println("  [2, 2]")
    println("-x: ", vps[1])
    println("+x: ", vps[2])
    println("-y: ", vps[3])
    println("+y: ", vps[4])

    println()
end


function test_veldt_setup_2d()
    print_test_header("Veldt Setup 2D")

    veldt1 = setup_veldt("./veldts/veldt1.yml")

    println("veldt1: ", veldt1)
end


function test_veldt_setup_3d()
    print_test_header("Veldt Setup 3D")

    veldt2 = setup_veldt("./veldts/veldt2.yml")

    println("veldt2: ", veldt2)
end


function test_simulation_2d()
    print_test_header("Veldt 2D Simulation")

    veldt1 = setup_veldt("./veldts/veldt3.yml")
    simulate(veldt1, 10)
end


function test_simulation_3d()
    print_test_header("Veldt 3D Simulation")

    veldt1 = setup_veldt("./veldts/veldt4.yml")
    simulate(veldt1, 10)
end


function test_diffusion_2d()
    print_test_header("Veldt 2D Diffusion")

    println("  * veldt5")
    veldt5 = setup_veldt("./veldts/veldt5.yml")
    simulate(veldt5, 20)

    println("  * veldt6")
    veldt6 = setup_veldt("./veldts/veldt6.yml")
    simulate(veldt6, 20)
end


function test_diffusion_3d()
    print_test_header("Veldt 3D Diffusion")

    println("  * veldt7")
    veldt7 = setup_veldt("./veldts/veldt7.yml")
    simulate(veldt7, 20)

    println("  * veldt8")
    veldt8 = setup_veldt("./veldts/veldt8.yml")
    simulate(veldt8, 20)
end


function test_simulation_setup_2d()
    print_test_header("Simulation Setup 2D")

    setup_simulation("./simulations/simulation1.yml")
end


function test_simulation_setup_3d()
    print_test_header("Simulation Setup 3D")

    setup_simulation("./simulations/simulation2.yml")
end


function test_substitution_matrix()
    print_test_header("SubstitutionMatrix")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Read SubstitutionMatrix")
    sm1 = read_substitution_matrix("./evolution_params/sub1.txt", chem1)

    println(sm1.substitutions['A'])
    println(sm1.substitutions['B'])
    println(sm1.substitutions['('])
    println(sm1.substitutions['/'])
end


function test_evolution_params()
    print_test_header("SelectionParams and MutationParams")

    sel_params, mut_params = read_evolution_params("./evolution_params/evolution_params1.yml")
end


function test_mutate()
    print_test_header("Mutate")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Read SubstitutionMatrix")
    sm1 = read_substitution_matrix("./evolution_params/sub1.txt", chem1)

    println("  * Setup MutationParams")
    mut_params = MutationParams(snv_rate=0.01,
                                insertion_rate=0.002,
                                deletion_rate=0.002,
                                duplication_rate=0.002,
                                inversion_rate=0.002,
                                translocation_rate=0.002,
                                crossing_over=true)

    println("  * Create Genomes")
    genomes = [Genome("genome$i", genome_string(500, chem1), chem1)
               for i in 1:5]

    for genome in genomes
        println(genome)
    end

    child_genomes = mutate(genomes, mut_params)

    for genome in child_genomes
        println(genome)
    end
end


function test_select_genomes()
    print_test_header("Select Genomes")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Setup SelectionParams")
    sel_params1 = SelectionParams(gene_count, select_count=4)
    sel_params2 = SelectionParams(genome_length, select_count=4)

    println("  * Create Genomes")
    genomes = [Genome("genome$i", genome_string(500 + rand(-5:5), chem1), chem1)
               for i in 1:8]

    for genome in genomes
        fitness = sel_params1.fitness_function(genome)
        @printf "%s: %d\n" genome.name fitness
    end

    selected_genomes = select_genomes(genomes, sel_params1)

    println("Selected - gene_count")
    for genome in selected_genomes
        println("  ", genome.name)
    end

    for genome in genomes
        fitness = sel_params2.fitness_function(genome)
        @printf "%s: %d\n" genome.name fitness
    end

    selected_genomes = select_genomes(genomes, sel_params2)

    println("Selected - genome_length")
    for genome in selected_genomes
        println("  ", genome.name)
    end
end


function test_select_cells()
    print_test_header("Select Cells")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Setup SelectionParams")
    sel_params = SelectionParams()

    println("  * Create Genomes")
    cells = [Cell(Genome("genome$i", genome_string(500, chem1), chem1))
             for i in 1:8]

    # for genome in genomes
    #     println(genome)
    # end

    selected_cells = select_cells(cells, sel_params)

    for cell in selected_cells
        println(cell)
    end
end


function test_genetic_algorithm1()
    print_test_header("Genetic Algorithm 1")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Setup SelectionParams")
    sel_params = SelectionParams(gene_count, select_count=10)

    println("  * Setup MutationParams")
    mut_params = MutationParams(snv_rate=0.01,
                                insertion_rate=0.0002,
                                deletion_rate=0.0002,
                                duplication_rate=0.0001,
                                inversion_rate=0.00004,
                                translocation_rate=0.00004,
                                crossing_over=false)

    println("  * Create Genomes")
    genomes = [Genome("genome$i", genome_string(500 + rand(-5:5), chem1), chem1)
               for i in 1:20]

    for i in 1:20
        println("  * Selection ", i)
        for genome in genomes
            fitness = sel_params.fitness_function(genome)
            @printf "%s: %d\n" genome.uuid fitness
        end

        selected_genomes = select_genomes(genomes, sel_params)

        println("Selected - gene_count")
        for genome in selected_genomes
            println("  ", genome.uuid)
        end

        println("  * Replication")
        new_genomes = [Genome(sel_genome.string, sel_genome.chemistry)
                       for sel_genome in selected_genomes]
        replicated_genomes = vcat(selected_genomes, new_genomes)

        println("  * Mutation")
        genomes = mutate(replicated_genomes, mut_params)
    end
end


function test_genetic_algorithm2()
    print_test_header("Genetic Algorithm 2")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Setup SelectionParams and MutationParams")
    sel_params, mut_params = read_evolution_params("./evolution_params/evolution_params1.yml")

    println("  * Create Genomes")
    genomes = [Genome("genome$i", genome_string(500 + rand(-5:5), chem1), chem1)
               for i in 1:20]

    for i in 1:20
        println("  * Selection ", i)
        for genome in genomes
            fitness = sel_params.fitness_function(genome)
            @printf "%s: %d\n" genome.uuid fitness
        end

        selected_genomes = select_genomes(genomes, sel_params)

        println("Selected - gene_count")
        for genome in selected_genomes
            println("  ", genome.uuid)
        end

        println("  * Replication")
        new_genomes = [Genome(sel_genome.string, sel_genome.chemistry)
                       for sel_genome in selected_genomes]
        replicated_genomes = vcat(selected_genomes, new_genomes)

        println("  * Mutation")
        genomes = mutate(replicated_genomes, mut_params)
    end
end

function test_metabolism1()
    print_test_header("Metabolism 1")

    println("  * Setup Chemistry")
    chem1 = setup_chemistry("./chemistries/chemistry1.yml")

    println("  * Reactions")
    reactions1 = [Reaction(["A", "B"], ["AB"], chem1),
                  Reaction(["B", "A"], ["BA"], chem1),
                  Reaction(["AB", "BA"], ["ABBA"], chem1),
                  Reaction(["ABBA"], ["ABB", "A"], chem1),
                  Reaction(["ABBA"], ["A", "BBA"], chem1)]
    println("  ** reactions1")
    for (i, reaction) in enumerate(reactions1)
        @printf "reaction%d: %s\n" i reaction
    end

    reactions2 = [Reaction(["A", "B"], ["AB"], chem1),
                  Reaction(["B", "A"], ["BA"], chem1),
                  Reaction(["A", "BB"], ["ABB"], chem1),
                  Reaction(["AB", "BA"], ["ABBA"], chem1),
                  Reaction(["AA"], ["A", "A"], chem1),
                  Reaction(["BB"], ["B", "B"], chem1),
                  Reaction(["ABBA"], ["ABB", "A"], chem1),
                  Reaction(["ABBA"], ["A", "BBA"], chem1)]
    println("  ** reactions2")
    for (i, reaction) in enumerate(reactions2)
        @printf "reaction%d: %s\n" i reaction
    end

    println("  * Metabolism")
    met1 = Metabolism(reactions1)
    met2 = Metabolism(reactions2)

    println("  ** met1")
    println("  *** nodes")
    for node in met1.nodes
        println(node)
    end
    @test length(met1.nodes) == 12

    println("  *** edges")
    for edge in met1.edges
        println(edge)
    end
    @test length(met1.edges) == 15

    submets1 = get_connected_submetabolisms(met1)
    println("  submetabolisms")
    @printf "submetabolisms: %d\n" length(submets1)
    @test length(submets1) == 1

    println("  ** met2")
    println("  *** nodes")
    for node in met2.nodes
        println(node)
    end
    @test length(met2.nodes) == 17

    println("  *** edges")
    for edge in met2.edges
        println(edge)
    end
    @test length(met2.edges) == 22

    submets2 = get_connected_submetabolisms(met1)
    println("  submetabolisms")
    @printf "submetabolisms: %d\n" length(submets2)
    @test length(submets2) == 1

    met3 = merge_metabolisms(met1, met2)
    println("  ** met3 = merge_metabolisms(met1, met2)")
    println("  *** nodes")
    for node in met3.nodes
        println(node)
    end
    @test length(met3.nodes) == 22

    println("  *** edges")
    for edge in met3.edges
        println(edge)
    end
    @test length(met3.edges) == 37

    submets3 = get_connected_submetabolisms(met3)
    println("  submetabolisms")
    @printf "submetabolisms: %d\n" length(submets3)
    @test length(submets3) == 1
end



function main()
    # Chemistry
    # test_element()
    # test_element_table()
    # test_bond()
    # test_bond_table()
    # test_chemistry()
    test_chemistry_setup()
    # test_molecule()
    # test_reaction()

    # Genome
    # test_gene()
    # test_genome()

    # Simulation
    # test_cell()
    # test_veldt_point()
    # test_veldt()
    # test_veldt_setup_2d()
    # test_veldt_setup_3d()
    # test_simulation_2d()
    # test_simulation_3d()
    # test_diffusion_2d()
    # test_diffusion_3d()
    # test_simulation_setup_2d()
    # test_simulation_setup_3d()

    # Evolution
    # test_substitution_matrix()
    # test_evolution_params()
    # test_mutate()
    # test_select_genomes()
    # test_select_cells()
    # test_genetic_algorithm1()
    # test_genetic_algorithm2()

    # Metabolism
    # test_metabolism1()
end

main()
