using MaxwellBase
using SimpleConstants
using Test
using Statistics: mean
using AbbreviatedTypes

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellBase" begin

include("enumtype.jl")
include("util.jl")
include("phys.jl")
include("field.jl")
include("material.jl")
include("object.jl")
include("smoothing.jl")
include("param.jl")

# end  # @testset "MaxwellBase"
