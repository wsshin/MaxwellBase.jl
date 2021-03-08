export create_field_array, field_arr2vec

# About the order of indices of fKd:
#
# Like paramKd in assignment.jl, we index fKd as fKd[i,j,k,w], where (i,j,k) are positional
# indices and w is the Cartesian component index.
#
# The reason for this choice is the same as paramKd.  For example, in assigning source
# values to jKd (an array of current density), we usually fix a component first and then
# assign the same value to a range of (i,j,k).  This can be more efficiently done with
# jKd[i,j,k,w], because a contiguous range of (i,j,k) actually corresponds to a contiguous
# memory block.
#
# I think I will need to index the E- and H-fields the same way for matrix-free operations.
# When we perform curl operation on these field arrays, we implement the ∂/∂w operation as
# fixing the component to differentiate first and then perform the differentiation.
# Therefore, again the E[i,j,k,w] indexing scheme results in an operation on a more
# contiguous block in memory space.

create_field_array(N::AbsVecInteger; ncmp::Int=3) = create_field_array(SInt{length(N)}(N), ncmp=ncmp)
create_field_array(N::SInt; ncmp::Int=3) = zeros(CFloat, N.data..., ncmp)


# Create a vector view of the field array.
#
# Note that reshape(fKd, :) creates a vector that shares the memory with fKd, i.e., if the
# entries of the resulting vector change, the corresponding entries of fKd change as well.
# Similarly, PermutedDimsArray(fKd, ...) creates a dimension-permuted array that shares the
# memory with fKd.  On the other hand, a related function permutedims(fKd, ...) creates a
# dimension-permuted array on a new memory space, so shouldn't be used.
field_arr2vec(fKd::AbsArrNumber{K₊₁};  # field array; K₊₁ = K+1, where K is space dimension and 1 is dimension for Cartesian components
              order_cmpfirst::Bool=true
              ) where {K₊₁} =
    order_cmpfirst ? reshape(PermutedDimsArray(fKd, (K₊₁, ntuple(identity,Val(K₊₁-1))...)), :) : reshape(fKd,:)
