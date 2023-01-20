module MaxwellBase

# @reexport makes all exported symbols of the exported packages available in module using
# MaxwellBase.  Basically, it is equivalent to calling `using ...` prepended by @reexport
# automatically upon calling `using MaxwellBase`.
using Reexport
@reexport using LinearAlgebra, SparseArrays, StaggeredGridCalculus, GeometryPrimitives
using ArgCheck
using Base.Threads: @threads  # used in param.jl
using AbbreviatedTypes; @define_types_with(Float64)
using SimpleConstants

export S²Number1, S²Number2, S²Number3, AbsMatParam, ParamInd, ObjInd

## Type aliases
const S²Number1 = S²Number{1,1}
const S²Number2 = S²Number{2,4}
const S²Number3 = S²Number{3,9}

const AbsMatParam = Union{Number,AbsVecNumber,AbsMatNumber}
const ParamInd = UInt8  # change this to handle more than 2⁸ = 256 materials
const ObjInd = UInt16  # change this to handle more than 2¹⁶ = 65536 objects


# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("enumtype.jl")
include("util.jl")
include("unit.jl")
include("material.jl")
include("object.jl")
include("field.jl")
include("assignment.jl")
include("smoothing.jl")
include("param.jl")

end # module MaxwellBase
