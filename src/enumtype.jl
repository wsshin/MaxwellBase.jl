export AbstractFieldType, Efield,Hfield

export gt4ft  # functions

abstract type AbstractFieldType <: AbstractBinaryType end
show(io::IO, ::MIME"text/plain", ft::AbstractFieldType) = show(io, ft)
broadcastable(ft::AbstractFieldType) = Ref(ft)

struct Efield <: AbstractFieldType end
show(io::IO, ::Efield) = print(io, "E-field")

getindex(A::Tuple, ::Efield) = first(A)
getindex(A::AbstractVector, ::Efield) = first(A)
setindex!(A, a, ::Efield) = setindex!(A, a, firstindex(A))  # (A[begin] = a; A)

conj(::Efield) = Hfield()

struct Hfield <: AbstractFieldType end
show(io::IO, ::Hfield) = print(io, "H-field")

getindex(A::Tuple, ::Hfield) = last(A)
getindex(A::AbstractVector, ::Hfield) = last(A)
setindex!(A, a, ::Hfield) = setindex!(A, a, lastindex(A))  # (A[end] = a; A)

conj(::Hfield) = Efield()

instances(::Type{AbstractFieldType}) = (Efield(), Hfield())


# Given boundary field types, determine whether the grid planes specified by the given field
# type ft are primal or dual grid planes in the Cartesian directions.
# Note that the boundary field type is the type of the fields tangential to the boundary.
# The fields defined on but normal to the boundary have the complementary type.
gt4ft(ft::AbstractFieldType, boundft::AbstractFieldType) = (ft==boundft ? PrimalGrid() : DualGrid())
