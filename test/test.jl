# John Eargle (mailto: jeargle at gmail.com)
# 2018
# test

using Uveldt


function test_element()
    println("***")
    println("*** Element")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el3 = Element('C', 3)

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

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()

    println("el_table1: ", el_table1)
    println()

    println("*** add Elements")
    add_elements(el_table1, [el1, el2])
    println("el_table1: ", el_table1)
    println("el_table1.elements[\"A\"]: ", el_table1.elements['A'])
    println("el_table1.elements[\"B\"]: ", el_table1.elements['B'])
    println()
end

function test_molecule()
    println("***")
    println("*** Molecule")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2])

    mol1 = Molecule("mol1", "AAA", el_table1)
    mol2 = Molecule("mol2", "BBB", el_table1)
    mol3 = Molecule("mol3", "ABA", el_table1)
    println("mol1: ", mol1)
    println("mol2: ", mol2)
    println("mol3: ", mol3)
    println()

    println("mass(mol1): ", mass(mol1))
    println("mass(mol2): ", mass(mol2))
    println("mass(mol3): ", mass(mol3))
    println()
end

function test_bond()
    println("***")
    println("*** Bond")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2])

    bond1 = Bond(el1, el1, el_table1, -5.5, 1.5)
    bond2 = Bond(el1, el2, el_table1, -4.5, 0.5)
    bond3 = Bond(el2, el2, el_table1, 3.5, 3.5)
    println("bond1: ", bond1)
    println("bond2: ", bond2)
    println("bond3: ", bond3)
    println()
end

function test_bond_table()
    println("***")
    println("*** BondTable")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2])

    bond1 = Bond(el1, el1, el_table1, -5.5, 1.5)
    bond2 = Bond(el1, el2, el_table1, -4.5, 0.5)
    bond3 = Bond(el2, el2, el_table1, 3.5, 3.5)

    b_table1 = BondTable(el_table1)
    println("b_table1: ", b_table1)

    println("*** add Bonds")
    add_bond(b_table1, bond1)
    add_bond(b_table1, bond2)
    add_bond(b_table1, bond3)
    println("b_table1: ", b_table1)

    println("*** get Bonds")
    println("get_bond(b_table1, 'A', 'A'): ", get_bond(b_table1, 'A', 'A'))
    println("get_bond(b_table1, 'B', 'B'): ", get_bond(b_table1, 'B', 'B'))
    println("get_bond(b_table1, 'A', 'B'): ", get_bond(b_table1, 'A', 'B'))
    println("get_bond(b_table1, 'B', 'A'): ", get_bond(b_table1, 'B', 'A'))
    println()
end

function test_gene()
    println("***")
    println("*** Gene")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el3 = Element('C', 3)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2, el3])

    genes = []
    push!(genes, Gene(1, "(A*B)", el_table1))
    push!(genes, Gene(21, "(A/B)", el_table1))
    push!(genes, Gene(45, "(A*B/C)", el_table1))
    push!(genes, Gene(132, "(A*B*C)", el_table1))
    push!(genes, Gene(154, "(A*//*A)", el_table1))
    push!(genes, Gene(154, "(A/*A*/A)", el_table1))
    push!(genes, Gene(154, "(*/AAA/*/*/)", el_table1))

    for i in 1:length(genes)
        @printf("gene%d: %s\n", i, genes[i])
    end

    println("*** Gene transcripts")
    for i in 1:length(genes)
        @printf("gene%d.transcript: %s\n", i, genes[i].transcript)
        reactants, products = parse_reaction(genes[i].transcript, el_table1)
        println("reactants: ", reactants)
        println("products: ", products)
    end

    println()
end

function test_genome()
    println("***")
    println("*** Genome")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2])

    genome_str = genome_string(200, el_table1)
    println("genome_string: ", genome_str)

    genome1 = Genome("genome1", genome_str, el_table1)
    println("genome1: ", genome1)

    genes = find_genes(genome1)

    println("*** find Genes")
    for gene in genes
        @printf("%d %s\n", gene.location, gene.string)
    end

    println("*** Gene transcripts")
    pseudogene_count = 0
    for i in 1:length(genes)
        if is_pseudogene(genes[i])
            println("*** pseudogene: " * genes[i].string)
            pseudogene_count += 1
        else
            @printf("gene%d.transcript: %s\n", i, genes[i].transcript)
            reactants, products = parse_reaction(genes[i].transcript, el_table1)
            println("reactants: ", reactants)
            println("products: ", products)
        end
    end

    @printf("%d good genes\n", length(genes)-pseudogene_count)
    @printf("%d pseudogenes\n", pseudogene_count)
    @printf("%d total\n", length(genes))

    println()
end


function test_reaction()
    println("***")
    println("*** Reaction")
    println("***")

    el1 = Element('A', 1)
    el2 = Element('B', 2)
    el_table1 = ElementTable()
    add_elements(el_table1, [el1, el2])

    bond1 = Bond(el1, el1, el_table1, -5.5, 1.5)
    bond2 = Bond(el1, el2, el_table1, -4.5, 0.5)
    bond3 = Bond(el2, el2, el_table1, 3.5, 3.5)

    b_table1 = BondTable(el_table1)

    println("*** add Bonds")
    add_bond(b_table1, bond1)
    add_bond(b_table1, bond2)
    add_bond(b_table1, bond3)
    println("b_table1: ", b_table1)

    println("*** Reactions")
    reaction1 = Reaction(["A", "A"], ["AA"], el_table1, b_table1)
    reaction2 = Reaction(["B", "B"], ["BB"], el_table1, b_table1)
    reaction3 = Reaction(["AB", "BA"], ["ABBA"], el_table1, b_table1)
    println("reaction1: ", reaction1)
    println("reaction2: ", reaction2)
    println("reaction3: ", reaction3)

    println()

end

function main()
    # test_element()
    # test_element_table()
    # test_molecule()
    # test_bond()
    test_bond_table()
    # test_gene()
    # test_genome()
    test_reaction()
end

main()
