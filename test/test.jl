# John Eargle (mailto: jeargle at gmail.com)
# 2018-2019
# test

using uveldt


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
    println("***")
    println("*** Element")
    println("***")

    el1, el2, el3 = create_elements(3)

    println("el1: ", el1)
    println("el2: ", el2)
    println("el3: ", el3)
    println()

    println("mass(el1): ", mass(el1))
    println("mass(el2): ", mass(el2))
    println("mass(el3): ", mass(el3))
    println()
end

function test_element_table()
    println("***")
    println("*** ElementTable")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    println("el_table1: ", el_table1)
    println("el_table1.elements[\"A\"]: ", el_table1.elements['A'])
    println("el_table1.elements[\"B\"]: ", el_table1.elements['B'])
    println()
end

function test_bond()
    println("***")
    println("*** Bond")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bond1, bond2, bond3 = create_bonds(el_table1, energies1)
    println("bond1: ", bond1)
    println("bond2: ", bond2)
    println("bond3: ", bond3)
    println()
end

function test_bond_table()
    println("***")
    println("*** BondTable")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)
    println("b_table1: ", b_table1)

    println("  * get Bonds")
    println("get_bond(b_table1, 'A', 'A'): ", get_bond(b_table1, 'A', 'A'))
    println("get_bond(b_table1, 'B', 'B'): ", get_bond(b_table1, 'B', 'B'))
    println("get_bond(b_table1, 'A', 'B'): ", get_bond(b_table1, 'A', 'B'))
    println("get_bond(b_table1, 'B', 'A'): ", get_bond(b_table1, 'B', 'A'))
    println()
end

function test_chemistry()
    println("***")
    println("*** Chemistry")
    println("***")

    elements1 = create_elements(2)
    el_table1 = ElementTable(elements1)

    elements2 = create_elements(3)
    el_table2 = ElementTable(elements2)

    println("  * add Bonds")
    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (-4.5, 0.5), (3.5, 3.5), (3.5, 3.5), (3.5, 3.5)]
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

end

function test_molecule()
    println("***")
    println("*** Molecule")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    mol1 = Molecule("mol1", "AAA", chem1)
    mol2 = Molecule("mol2", "BBB", chem1)
    mol3 = Molecule("mol3", "ABA", chem1)
    println("mol1: ", mol1)
    println("mol2: ", mol2)
    println("mol3: ", mol3)
    println()

    println("mass(mol1): ", mass(mol1))
    println("mass(mol2): ", mass(mol2))
    println("mass(mol3): ", mass(mol3))
    println()
end

function test_reaction()
    println("***")
    println("*** Reaction")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
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
        @printf("reaction%d: %s\n", i, reaction)
    end

    println("  * Pores")
    reactions = [Reaction("A", chem1),
                 Reaction("BB", chem1),
                 Reaction("ABBA", chem1)]
    for (i, reaction) in enumerate(reactions)
        @printf("pore%d: %s\n", i, reaction)
    end

    println("  * Transporters")
    reactions = [Reaction("A", chem1, transport_in, -3.0),
                 Reaction("BB", chem1, transport_in, -4.0),
                 Reaction("ABBA", chem1, transport_in, -5.0),
                 Reaction("A", chem1, transport_out, -3.0),
                 Reaction("BB", chem1, transport_out, -4.0),
                 Reaction("ABBA", chem1, transport_out, -5.0)]
    for (i, reaction) in enumerate(reactions)
        @printf("pore%d: %s\n", i, reaction)
    end

    println()

end

function test_gene()
    println("***")
    println("*** Gene")
    println("***")

    elements = create_elements(3)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (-4.5, 0.5), (3.5, 3.5), (3.5, 3.5), (3.5, 3.5)]
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
        @printf("gene%d: %s\n", i, genes[i])
    end

    println("  * Gene transcripts")
    reactions = []
    for i in 1:length(genes)
        @printf("gene%d.transcript: %s\n", i, genes[i].transcript)
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
        @printf("reaction%d: %s\n", i, reactions[i])
    end

    println()
end

function test_genome()
    println("***")
    println("*** Genome")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
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
    genome_str2 = genome_info[1][2]
    println("genome_str2: ", genome_str2)
    println("name1 == name2: ", genome_info[1][1] == "genome1")
    println("genome_str1 == genome_str2: ", genome_str1 == genome_str2)

    genes = find_genes(genome1)

    println("  * find Genes")
    for gene in genes
        @printf("%d %s\n", gene.location, gene.string)
    end

    println("  * Gene transcripts")
    pseudogene_count = 0
    reaction_count = 0
    pore_count = 0
    transport_in_count = 0
    transport_out_count = 0
    for i in 1:length(genes)
        @printf("gene%d.transcript: %s\n", i, genes[i].transcript)
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
            elseif reaction_type == transport_in
                transport_in_count += 1
            elseif reaction_type == transport_out
                transport_out_count += 1
            end
        end
    end

    @printf("%d good genes\n", length(genes)-pseudogene_count)
    @printf("  %d reactions\n", reaction_count)
    @printf("  %d pores\n", pore_count)
    @printf("  %d transport in\n", transport_in_count)
    @printf("  %d transport out\n", transport_out_count)
    @printf("%d pseudogenes\n", pseudogene_count)
    @printf("%d total\n", length(genes))

    add_snps!(genome1, 0.01)
    add_insertions!(genome1, 0.002)
    remove_deletions!(genome1, 0.002)

    println()
end

function test_cell()
    println("***")
    println("*** Cell")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    genome_str = genome_string(200, chem1)
    println("genome_string: ", genome_str)

    genome1 = Genome("genome1", genome_str, chem1)
    println("genome1: ", genome1)

    println("  * Cells")
    cell1 = Cell(genome1)
    println("cell1: ", cell1)

    println("  * add Molecule counts")
    cell1.molecule_counts[1]["AAA"] = 33
    cell1.molecule_counts[1]["BBB"] = 44
    println("cell1: ", cell1)

    println()
end

function test_veldt_point()
    println("***")
    println("*** VeldtPoint")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    genome_str = genome_string(200, chem1)
    genome1 = Genome("genome1", genome_str, chem1)
    cell1 = Cell(genome1)
    cell1.molecule_counts[1]["AAA"] = 33
    cell1.molecule_counts[1]["BBB"] = 44

    println("  * VeldtPoints")
    veldt_pt1 = VeldtPoint()
    veldt_pt2 = VeldtPoint(molecules=["AAA", "BBB", "ABA"])
    veldt_pt3 = VeldtPoint(cell=cell1)
    println("veldt_pt1: ", veldt_pt1)
    println("veldt_pt2: ", veldt_pt2)
    println("veldt_pt3: ", veldt_pt3)

    println("  * add Molecule counts")
    veldt_pt1.molecule_counts[1]["AAA"] = 300
    veldt_pt1.molecule_counts[1]["BBB"] = 400
    veldt_pt2.molecule_counts[1]["AAA"] = 330
    veldt_pt2.molecule_counts[1]["BBB"] = 440
    veldt_pt2.molecule_counts[1]["AAA"] = 333
    veldt_pt2.molecule_counts[1]["BBB"] = 444
    println("veldt_pt1: ", veldt_pt1)
    println("veldt_pt2: ", veldt_pt2)
    println("veldt_pt3: ", veldt_pt3)

    println()
end

function test_veldt()
    println("***")
    println("*** Veldt")
    println("***")

    elements = create_elements(2)
    el_table1 = ElementTable(elements)

    energies1 = [(-5.5, 1.5), (-4.5, 0.5), (3.5, 3.5)]
    bonds = create_bonds(el_table1, energies1)

    b_table1 = BondTable(bonds)

    println("  * create Chemistry")
    chem1 = Chemistry(el_table1, b_table1)

    genome_str1 = genome_string(200, chem1)
    genome1 = Genome("genome1", genome_str1, chem1)
    genome_str2 = genome_string(200, chem1)
    genome2 = Genome("genome2", genome_str2, chem1)

    cell1 = Cell(genome1)
    cell1.molecule_counts[1]["AAA"] = 11
    cell1.molecule_counts[1]["BBB"] = 22
    cell2 = Cell(genome2)
    cell2.molecule_counts[1]["AA"] = 33
    cell2.molecule_counts[1]["AB"] = 44

    println("  * Veldt")
    veldt1 = Veldt([3, 4])
    println("veldt1: ", veldt1)

    println("  * initialize Molecule counts")
    init_molecules(veldt1, [1, 1], Dict("AAA" => 10, "BBB" => 10))
    init_molecules(veldt1, [1, 2], Dict("AAA" => 10, "BBB" => 20))
    init_molecules(veldt1, [1, 3], Dict("AAA" => 10, "BBB" => 30))
    init_molecules(veldt1, [3, 3], Dict("AAA" => 30, "BBB" => 30))
    println("veldt1: ", veldt1)

    println("  * add Cells")
    add_cell(veldt1, [1, 1], cell1)
    add_cell(veldt1, [2, 2], cell2)
    println("veldt1: ", veldt1)

    println()
end



function main()
    # test_element()
    # test_element_table()
    # test_bond()
    # test_bond_table()
    # test_chemistry()
    # test_molecule()
    # test_reaction()
    # test_gene()
    test_genome()
    # test_cell()
    # test_veldt_point()
    # test_veldt()
end

main()
