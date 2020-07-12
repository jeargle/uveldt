# Run create_sysimage() from this directory.

using Distributions
using Printf
using Random
using YAML

str1 = "test string"
location1 = rand(0:length(str1))
str2 = randstring("abcdefg", 30)
rate = 0.01
geom_dist = Geometric(rate)
location2 = 1 + rand(geom_dist)


println()
@printf "PackageCompiler so Builder\n"
@printf "  str1: %s\n" str1
@printf "  str2: %s\n" str2
@printf "  location1: %d\n" location1
@printf "  location2: %d\n" location2

setup = YAML.load(open("chemistries/chemistry1.yml"))
if haskey(setup, "elements")
    for el_info in setup["elements"]
        name = el_info["name"][1]
        mass = el_info["mass"]
        @printf "  element: %s %d\n" name mass
    end
end
