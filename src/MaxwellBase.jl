module MaxwellBase

# @reexport makes all exported symbols of the exported packages available in module using
# MaxwellBase.  Basically, it is equivalent to calling `using ...` prepended by @reexport
# automatically upon calling `using MaxwellBase`.
using Reexport
@reexport using LinearAlgebra, SparseArrays, StaggeredGridCalculus, GeometryPrimitives
using StaticArrays
using Base.Threads: @threads  # used in param.jl
using SimpleConstants
using ArgCheck

import Base: show, broadcastable, getindex, setindex!, conj, instances

export AbstractMatParam, ParamInd, ObjInd

## Type aliases
const AbstractMatParam = Union{Number,AbstractVector{<:Number},AbstractMatrix{<:Number}}
const ParamInd = UInt8  # change this to handle more than 2⁸ = 256 materials
const ObjInd = UInt16  # change this to handle more than 2¹⁶ = 65536 objects


# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("binarytype.jl")
include("util.jl")
include("unit.jl")
include("material.jl")
include("object.jl")
include("field.jl")
include("assignment.jl")
include("smoothing.jl")
include("param.jl")

end # module MaxwellBase
