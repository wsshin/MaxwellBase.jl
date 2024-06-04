using MaxwellBase
using StaticArrays
using SimpleConstants
using Test
using Statistics: mean

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

# @testset "MaxwellBase" begin

include("enumtype.jl")
include("util.jl")
include("field.jl")
include("material.jl")
include("object.jl")
include("smoothing.jl")
include("param.jl")

# end  # @testset "MaxwellBase"
