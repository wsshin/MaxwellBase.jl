export FieldType, BC  # types
export EH, nE, nH, nEH  # instances
export ft2gt  # functions

# Field types
const nE, nH = 1, 2  # E-, H-fields
const nEH = SVec(nE, nH)
@enum FieldType EE=nE HH
const EH = SVec(EE, HH)
for ins in instances(FieldType); @eval export $(Symbol(ins)); end  # export all instances
Base.string(ins::FieldType) = ins==EE ? "E" : "H"
StaggeredGridCalculus.alter(ins::FieldType) = ins==EE ? HH : EE

# Given boundary field types, determine whether the grid planes specified by the given field
# type ft are primal or dual grid planes in the Cartesian directions.
# Note that the boundary field type is the type of the fields tangential to the boundary.
# The fields defined on but normal to the boundary have the complementary type.
ft2gt(ft::FieldType, boundft::FieldType) = PD[2 - (boundft==ft)]

# Boundary conditions
@enum BC PERIODIC=1 CONDUCTING  # periodic boundary condition, perfect conductor boundary condition
for ins in instances(BC); @eval export $(Symbol(ins)); end  # export all instances
Base.string(ins::BC) = ins==PERIODIC ? "periodic" : "perfect conductor"
