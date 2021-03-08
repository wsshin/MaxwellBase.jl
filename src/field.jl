export create_field_array, field_arr2vec, isfield_ortho_shape


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


# Determine if the field subspace is orthogonal to the shape subspace; used in smoothing.jl
# and source/*.jl.
#
# The shape and field subspaces are the subspaces of the 3D space where the shapes and
# fields lie.  In standard 3D problems, both subspaces are 3D.  However there are other
# cases as well: for example, in 2D TM problems, the shapes are on the xy-plane (2D space),
# but the E-fields are along the z-direction (1D space).
#
# For the shape and field subspace dimensions K and Kf, orthogonal subspaces can occur only
# when K + Kf ≤ 3, because if two subspaces are orthogonal to each other and K + Kf > 3, we
# can choose K + Kf linearly independent vectors in the direct sum of the two subspaces,
# which is impossible as the direct sum should still be a subspace of the 3D space.
#
# Considering the above, there are only a few combinations of K and Kf that allow orthogonal
# subspaces: (K,Kf) = (2,1), (1,1), (1,2).
# - (K,Kf) = (2,1).  This happens in 2D TE or TM problems.  The shape subspace is the 2D xy-
# plane, but the magnetic (electric) field subspace in TE (TM) problems is the 1D z-axis.
# - (K,Kf) = (1,1).  This happens in the 1D TEM problems with isotropic materials.  The
# shape subspace is the 1D z-axis, but the E- and H-field spaces are the 1D x- and y-axes.
# - (K,Kf) = (1,2).  This happens in the 1D TEM problem with anisotropic materials.
#
# Note that we can always solve problems as if they are 3D problems.  So, the consideration
# of the cases with K, Kf ≠ 3 occurs only when we can build special equations with reduced
# number of degrees of freedom, like in the TE, TM, TEM equations.  In such cases, we find
# that the two subspaces are always orthogonal if K ≠ Kf.  In fact, smooth_param!() is
# written such that the Kottke's subpixel smoothing algorithm that decomposes the field into
# the components tangential and normal to the shape surface is applied only when the shape
# and field subspaces coincide.  (This makes sense, because the inner product between the
# field and the direction normal that needs to be performed to decompose the field into such
# tangential and normal components is defined only when the field and the direction normal
# are in the same vector space.)  Therefore, if we want to apply Kottke's subpixel smoothing
# algorithm that depends on the surface normal direction, we have to construct the problem
# such that K = Kf (but the converse is not true: K = Kf does not imply the use of Kottke's
# subpixel smoothing algorithm that depends on the surface normal direction; in other words,
# the case with K = Kf can still use the subpixel smoothing algorithm that assumes the
# orthogonality between the shape and field subspace, such that the field is always
# tangential to the shape surface.)
#
# The contraposition of the above statement is that if K ≠ Kf, then Kottke's subpixel
# smoothing algorithm that does NOT depend on the direction normal is applied.  This
# subpixel smoothing algorithm assumes that the field subspace is orthogonal to the shape
# subspace, such that the field has only the tangential component to the surface.
# Therefore, if we pass K ≠ Kf to smooth_param!(), we should make sure that the field
# subspace is orthogonal to the shape subspace.  Note that this does not exclude the
# possibility of K ≠ Kf while the two subspaces are nonorthogonal; it just means that we
# must formulate the problem differently in such cases by, e.g., decomposing the equations,
# in order to use smooth_param!().
#
# As noted earlier, K = Kf can stil include cases where the shape and field subspaces are
# orthogonal.  Because K + Kf = 2K should be ≤ 3 in+ such cases, we conclude that K = Kf = 1
# is the only case where the two subspaces could be orthogonal while K = Kf.  In fact, when
# K = Kf, we will assume that the two subspaces are orthogonal, because the problems with 1D
# slabs and the field along the slab thickness direction are not interesting from the EM
# wave propagation perspectives.
isfield_ortho_shape(Kf, K) = Kf≠K || Kf==1
