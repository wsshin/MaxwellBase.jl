module MaxwellBase

# @reexport makes all exported symbols of the exported packages available in module using
# MaxwellBase.  Basically, it is equivalent to calling `using ...` prepended by @reexport
# automatically upon calling `using MaxwellBase`.
using Reexport
@reexport using LinearAlgebra, SparseArrays, StaggeredGridCalculus, GeometryPrimitives
using AbbreviatedTypes
using Base.Threads: @threads  # used in param.jl

export SSComplex1, SSComplex2, SSComplex3, ParamInd, ObjInd

## Type aliases
# Below, use Int instead of Int64 for compatibility with 32-bit systems (e.g., x86 in appveyor.yml).
const SSComplex{K,L} = SMatrix{K,K,CFloat,L}

const SSComplex1 = SSComplex{1,1}
const SSComplex2 = SSComplex{2,4}
const SSComplex3 = SSComplex{3,9}

const MatParam = Union{Number,AbsVecNumber,AbsMatNumber}
const ParamInd = UInt8  # change this to handle more than 2⁸ = 256 materials
const ObjInd = UInt16  # change this to handle more than 2¹⁶ = 65536 objects


# The order of inclusion matters: if types or functions in file A are used in file B, file A
# must be included first.
include("enumtype.jl")
include("util.jl")
include("phys.jl")
include("material.jl")
include("object.jl")
include("field.jl")
include("assignment.jl")
include("smoothing.jl")
include("param.jl")

end # module MaxwellBase
