# John Eargle (mailto: jeargle at gmail.com)
# 2018
# test

using Uveldt


function test_element()
    println("***")
    println("*** Element")
    println("***")

    el1 = Element("A", 1)
    el2 = Element("B", 2)
    el3 = Element("C", 3)

    println(el1)
    println(el2)
    println(el3)
    println()
end

function test_element_table()
    println("***")
    println("*** ElementTable")
    println("***")

    el1 = Element("A", 1)
    el2 = Element("B", 2)
    el_table1 = ElementTable()

    println("el_table1: ", el_table1)
    println()

    println("*** add Elements")
    add_element(el_table1, el1)
    add_element(el_table1, el2)
    println("el_table1: ", el_table1)
    println("el_table1.elements[\"A\"]: ", el_table1.elements["A"])
    println("el_table1.elements[\"B\"]: ", el_table1.elements["B"])
    println()
end



function main()
    test_element()
    test_element_table()
end

main()
