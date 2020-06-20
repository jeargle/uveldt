using Distributions
using Printf
using Random

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
