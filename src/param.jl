export create_paramop, apply_paramop!

create_paramop(paramKd::AbstractArray{<:Number},  # material parameter array
               isfwd_in::AbstractVector{Bool},  # isfwd_in[w] = true if w-component of input field is forward-averaged to be interpolated at voxel corners
               ∆l::NTuple{K,AbstractVector{<:Number}},  # line segments to multiply with (defined at field voxel edges); vectors of length N
               ∆l′⁻¹::NTuple{K,AbstractVector{<:Number}},  # inverse of line segments to divide by (defined at field voxel corners); vectors of length N
               isbloch::AbstractVector{Bool},  # for K = 3, boundary conditions in x, y, z
               e⁻ⁱᵏᴸ::AbstractVector{<:Number}=ones(Int,K);  # for K = 3, Bloch phase factor in x, y, z
               cmpfirst::Bool=true  # true to use Cartesian-component-major ordering for more tightly banded matrix
               ) where {K} =
    create_paramop(paramKd, SVector{K}(isfwd_in), ∆l, ∆l′⁻¹, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ); cmpfirst)

# Version for scalar material parameters
function create_paramop(paramKd::AbstractArray{<:Number,K₊₂};  # scalar material parameter array
                        cmpfirst::Bool=true  # true to use Cartesian-component-major ordering for more tightly banded matrix
                        ) where {K₊₂}
    K = K₊₂ - 2
    N = SVector(ntuple(k->size(paramKd,k)-1, Val(K)))  # for -1, see comments about size of paramKd in assignment.jl
    kdiag = 0
    paramop = create_paramop_supdiag(paramKd, kdiag, N; cmpfirst)

    return paramop
end

# Below, N can be retrieved from the size of paramKd or ∆l, but I pass it to make the size
# more explicit.  This also makes the argument list similar to create_mean()'s, which
# includes N because it allows omitting ∆l and ∆l′ and using ∆l = ∆l′ = ones.((N...,)) as
# the default arguments.
function create_paramop(paramKd::AbstractArray{<:Number,K₊₂},  # material parameter array
                        isfwd_in::SVector{K,Bool},  # isfwd_in[w] = true if w-component of input field is forward-averaged to be interpolated at voxel corners
                        ∆l::NTuple{K,AbstractVector{<:Number}},  # line segments to multiply with (defined at field voxel edges); vectors of length N
                        ∆l′⁻¹::NTuple{K,AbstractVector{<:Number}},  # inverse of line segments to divide by (defined at field voxel corners); vectors of length N
                        isbloch::SVector{K,Bool},  # for K = 3, boundary conditions in x, y, z
                        e⁻ⁱᵏᴸ::SVector{K,<:Number};  # for K = 3, Bloch phase factor in x, y, z
                        cmpfirst::Bool=true  # true to use Cartesian-component-major ordering for more tightly banded matrix
                        ) where {K,K₊₂}
    @assert K₊₂==K+2

    Kf = size(paramKd, K+1)  # field dimension
    @assert size(paramKd,K+2)==Kf

    N = SVector(length.(∆l))
    @argcheck length.(∆l′⁻¹)==N.data  # corresponding entries of ∆l and ∆l′⁻¹ should have same length
    M = prod(N)  # number voxels

    # Make sure Kf ≥ 2; the material parameters should be tensorial.
    @argcheck Kf≥2  # paramKd should represent array of tensorial material parameters

    # Following Oskooi et al.'s 2009 Optics Letters paper, off-diagonal entries of
    # material parameter tensors (e.g., ε) are evaluated at the corners of voxels whose
    # edges are field lines (e.g., E).  Therefore, if w-normal voxel faces are primal
    # grid planes, the w-component of the input fields need to be averaged in the
    # backward direction to be interpolated at voxel corners.
    #
    # Once these interpolated input fields (e.g., Ex) are multiplied with off-diagonal
    # entries (e.g., εzx) of the material parameter tensor, we get the output fields
    # (e.g., Dzx of Dz = Dzx + Dzy + Dzz = εzx Ex + εzy Ey + εzz Ez) in the direction
    # normal (say v-direction) to the input fields.  These output fields are still at
    # the voxel corners.  Now, if the v-normal voxel faces are primal grid planes, the
    # v-component of the output fields need to be averaged in the forward direction to
    # be interpolated at voxel edges.
    #
    # In summary,
    #
    # - to obtain the fields to feed to the off-diagonal entries of material parameter
    # tensors, the w-component of the input fields need to be forward(backward)-averaged
    # along the w-direction if the w-normal voxel faces are dual (primal) grid planes,
    #
    # - to distribute the resulting output fields in the v(≠w)-direction back to voxel
    # edges, the output fields need to be forward(backward)-averaged along the
    # v-direction if the v-normal voxel faces are primal (dual) grid planes.
    #
    # This means if the w-boundaries are E-field boundaries,
    # - the input E(H)w-field needs to be backward(forward)-averaged in the w-direction,
    # - the output E(H)w-field needs to be forward(backward)-averaged in the w-direction.
    isfwd_out = .!isfwd_in  # SVector{K,Bool}; true if output fields need to be backward-averaged

    # For the output averaging, ∆l and ∆l′ are not supplied to create_mean in order to
    # create a simple arithmetic averaging operator.  This is because the area factor
    # matrix multiplied for symmetry to the left of the material parameter matrix
    # multiplies the same area factor to the two fields being averaged.  (See my notes
    # on Jul/16/2018 in MaxwellFDFD in Agenda.)
    Mout = create_mean(isfwd_out, N, isbloch, e⁻ⁱᵏᴸ; cmpfirst)

    # For the input averaging, ∆l and ∆l′ are supplied to create_mean in order to create
    # a line integral averaging operator.  This is because the inverse of the length
    # factor matrix multiplied for symmetry to the right of the material parameter
    # matrix divides the two fields being averaged by different (= nonuniform) line
    # segments.  The ∆l factors multiplied inside create_minfo cancel the effect of this
    # multiplication with the nonuniform line segments.  (See my notes on Jul/16/2018 in
    # MaxwellFDFD in Agenda.)
    Min = create_mean(isfwd_in, ∆l, ∆l′⁻¹, isbloch, e⁻ⁱᵏᴸ; cmpfirst)

    # Off-diagonal components of material parameter tensor
    # Min and Mout are block-diagonal matrices, but create_paramop_supdiag() generates
    # block-off-diagonal matrices for 1 ≤ kdiag ≤ Kf-1.  See my bullet point entitled
    # [Update (May/13/2018)] in RN - Subpixel Smoothing.
    kdiag = 1
    paramop = create_paramop_supdiag(paramKd, kdiag, N; cmpfirst)
    for kdiag = 2:Kf-1  # index of diagonal of material parameter tensor
        paramop += create_paramop_supdiag(paramKd, kdiag, N; cmpfirst)
    end
    paramop = Mout * paramop * Min

    # Diagonal components of material parameter tensor
    kdiag = 0
    paramop += create_paramop_supdiag(paramKd, kdiag, N; cmpfirst)  # diagonal components of material parameter tensor

    return paramop
end

# Create the kth superdiagonal components of the parameter operator (0 ≤ k ≤ Kf-1)
#
# About the meaning of kdiag
#
# Suppose material parameter tensor is Kf-by-Kf.  Different values of kdiag take different
# entries of the material parameter tensor and used to create a matrix.  For Kf = 3, here
# are the entries (indicated by X) taken for the three different values of kdiag.
#
# kdiag = 0:
# ⎡X O O⎤
# ⎢O X O⎥
# ⎣O O X⎦
#
# kdiag = 1:
# ⎡O X O⎤
# ⎢O O X⎥
# ⎣X O O⎦
#
# kdiag = 2 (= Kf-1):
# ⎡O O X⎤
# ⎢X O O⎥
# ⎣O X O⎦
#
# So, kdiag is the index counted from the main diagonal (kdiag = 0) along the superdiagonal
# direction.  See RN - Subpixel Smoothing > [Update (May/13/2018)].  Note that if the X is
# the vw-block, then it is constructed with the vw-component of the material parameter (e.g.,
# for ε, if X is the xy-block, then it is constructed with εxy).  Therefore, whether the
# direction of the superdiagonal shift is along the row- or column-direction is not relevant.
function create_paramop_supdiag(paramKd::AbstractArray{P,K₊₂},  # size of last two dimensions is Kf-by-Kf
                                kdiag::Integer,  # 0 ≤ kdiag ≤ Kf-1; index of diagonal of material parameter tensor to set (kdiag = 0: main diagonal)
                                N::SVector{K,Int};  # size of grid
                                cmpfirst::Bool=true  # true to use Cartesian-component-major ordering for more tightly banded matrix
                                ) where {P,K,K₊₂}
    @assert K₊₂==K+2

    Kf = size(paramKd, K+1)  # field dimension
    @assert size(paramKd,K+2)==Kf

    # Note that paramKd's i, j, k indices run from 1 to N+1 rather than to N, so we should
    # not iterate those indices from 1 to end (= N+1).
    M = prod(N)  # number of voxels
    KfM = Kf * M
    I = Vector{Int}(undef, KfM)
    J = Vector{Int}(undef, KfM)
    V = Vector{P}(undef, KfM)
    n = 0

    CI = CartesianIndices(N.data)
    LI = LinearIndices(N.data)
    for nv = 1:Kf  # row index of material parameter tensor
        istr, ioff = cmpfirst ? (Kf, nv-Kf) : (1, M*(nv-1))  # (row stride, row offset)
        nw = mod1(nv+kdiag, Kf)  # column index of material parameter tensor
        jstr, joff = cmpfirst ? (Kf, nw-Kf) : (1, M*(nw-1))  # (column stride, column offset)
        for ci = CI
            n += 1
            @inbounds ind = LI[ci]  # linear index of Yee's cell

            I[n] = istr * ind + ioff
            J[n] = jstr * ind + joff
            @inbounds V[n] = paramKd[ci,nv,nw]
        end
    end

    return sparse(I, J, V, KfM, KfM)
end

apply_paramop!(G::AbstractArray{<:Number},  # output field; in 3D, G[i,j,k,w] is w-component of G at (i,j,k)
               F::AbstractArray{<:Number},  # input field; in 3D, F[i,j,k,w] is w-component of F at (i,j,k)
               G′::AbstractArray{<:Number},  # first temporary storage for intermediate output; in 3D, G′[i,j,k,w] is w-component of G′ at (i,j,k)
               G′′::AbstractArray{<:Number},  # second temporary storage for intermediate output; in 3D, G′[i,j,k,w] is w-component of G′ at (i,j,k)
               paramKd::AbstractArray{<:Number},  # material parameter array
               isfwd_in::AbstractVector{Bool},  # isfwd_in[w] = true if w-component of input field is forward-averaged to be interpolated at voxel corners
               ∆l::NTuple{K,AbstractVector{<:Number}},  # line segments to multiply with (defined at field voxel edges); vectors of length N
               ∆l′⁻¹::NTuple{K,AbstractVector{<:Number}},  # line segments to divide by (defined at field voxel corners); vectors of length N
               isbloch::AbstractVector{Bool},  # for K = 3, boundary conditions in x, y, z
               e⁻ⁱᵏᴸ::AbstractVector{<:Number}=ones(1,K);  # for K = 3, Bloch phase factor in x, y, z
               α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
               overwrite::Bool=true  # false to add result to output storage G
               ) where {K} =
apply_paramop!(G, F, G′, G′′, paramKd, SVector{K}(isfwd_in), ∆l, ∆l′⁻¹, SVector{K}(isbloch), SVector{K}(e⁻ⁱᵏᴸ); α, overwrite)


# Version for scalar material parameters
function apply_paramop!(G::AbstractArray{<:Number,K₊₁},  # output field; in 3D, G[i,j,k,w] is w-component of G at (i,j,k)
                        F::AbstractArray{<:Number,K₊₁},  # input field; in 3D, F[i,j,k,w] is w-component of F at (i,j,k)
                        paramKd::AbstractArray{<:Number,K₊₂};  # material parameter array
                        α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
                        overwrite::Bool=true  # false to add result to output storage G
                        ) where {K₊₁,K₊₂}
    kdiag = 0
    apply_paramop_supdiag!(G, F, paramKd, kdiag; α, overwrite)  # diagonal components of material parameter tensor

    return nothing
end

function apply_paramop!(G::AbstractArray{<:Number,K₊₁},  # output field; in 3D, G[i,j,k,w] is w-component of G at (i,j,k)
                        F::AbstractArray{<:Number,K₊₁},  # input field; in 3D, F[i,j,k,w] is w-component of F at (i,j,k)
                        G′::AbstractArray{<:Number,K₊₁},  # first temporary storage for intermediate output; in 3D, G′[i,j,k,w] is w-component of G′ at (i,j,k)
                        G′′::AbstractArray{<:Number,K₊₁},  # second temporary storage for intermediate output; in 3D, G′[i,j,k,w] is w-component of G′ at (i,j,k)
                        paramKd::AbstractArray{<:Number,K₊₂},  # material parameter array
                        isfwd_in::SVector{K,Bool},  # isfwd_in[w] = true if w-component of input field is forward-averaged to be interpolated at voxel corners
                        ∆l::NTuple{K,AbstractVector{<:Number}},  # line segments to multiply with (defined at field voxel edges); vectors of length N
                        ∆l′⁻¹::NTuple{K,AbstractVector{<:Number}},  # line segments to divide by (defined at field voxel corners); vectors of length N
                        isbloch::SVector{K,Bool},  # for K = 3, boundary conditions in x, y, z
                        e⁻ⁱᵏᴸ::SVector{K,<:Number};  # for K = 3, Bloch phase factor in x, y, z
                        α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
                        overwrite::Bool=true  # false to add result to output storage G
                        ) where {K,K₊₁,K₊₂}
    @assert K₊₁==K+1 && K₊₂==K+2

    Kf = size(paramKd, K+1)  # field dimension
    @assert size(paramKd,K+2)==Kf

    # Make sure Kf ≥ 2; the material parameters should be tensorial.
    @argcheck Kf≥2  # paramKd should represent array of tensorial material parameters

    isfwd_out = .!isfwd_in  # SVector{K,Bool}; true if output fields need to be backward-averaged

    apply_mean!(G′, F, isfwd_in, ∆l, ∆l′⁻¹, isbloch, e⁻ⁱᵏᴸ)

    # Off-diagonal components of material parameter tensor
    kdiag = 1
    apply_paramop_supdiag!(G′′, G′, paramKd, kdiag)
    for kdiag = 2:Kf-1  # index of diagonal of material parameter tensor
        apply_paramop_supdiag!(G′′, G′, paramKd, kdiag, overwrite=false)
    end
    apply_mean!(G, G′′, isfwd_out, isbloch, e⁻ⁱᵏᴸ; α, overwrite)

    # Diagonal components of material parameter tensor
    kdiag = 0
    apply_paramop_supdiag!(G, F, paramKd, kdiag; α, overwrite=false)

    return nothing
end

# 3D
function apply_paramop_supdiag!(G::AbstractArray{<:Number,4},  # output field; G[i,j,k,w] is w-component of G at (i,j,k)
                                F::AbstractArray{<:Number,4},  # input field; F[i,j,k,w] is w-component of F at (i,j,k)
                                paramKd::AbstractArray{<:Number,5},  # size of last two dimensions is Kf-by-Kf
                                kdiag::Integer;  # 0 ≤ kdiag ≤ Kf-1; index of diagonal of material parameter tensor to set (kdiag = 0: main diagonal)
                                α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
                                overwrite::Bool=true)  # false to add result to output storage G
    Nx, Ny, Nz, Kf = size(G)  # Kf: field dimension
    @assert size(paramKd)[end] == size(paramKd)[end-1] == Kf

    overwrite && (G .= 0)
    for nv = 1:Kf  # row index of material parameter tensor
        Gᵥ = @view G[:,:,:,nv]
        nw = mod1(nv+kdiag, Kf)  # column index of material parameter tensor
        @threads for k = 1:Nz
            for j = 1:Ny, i = 1:Nx
                @inbounds Gᵥ[i,j,k] += paramKd[i,j,k,nv,nw] * F[i,j,k,nw]
            end
        end
    end

    return nothing
end

# 2D
function apply_paramop_supdiag!(G::AbstractArray{<:Number,3},  # output field; G[i,j,w] is w-component of G at (i,j)
                                F::AbstractArray{<:Number,3},  # input field; F[i,j,w] is w-component of F at (i,j)
                                paramKd::AbstractArray{<:Number,4},  # size of last two dimensions is Kf-by-Kf
                                kdiag::Integer;  # 0 ≤ kdiag ≤ Kf-1; index of diagonal of material parameter tensor to set (kdiag = 0: main diagonal)
                                α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
                                overwrite::Bool=true)  # false to add result to output storage G
    Nx, Ny, Kf = size(G)  # Kf: field dimension
    @assert size(paramKd)[end] == size(paramKd)[end-1] == Kf

    overwrite && (G .= 0)
    for nv = 1:Kf  # row index of material parameter tensor
        Gᵥ = @view G[:,:,nv]
        nw = mod1(nv+kdiag, Kf)  # column index of material parameter tensor
        @threads for j = 1:Ny
            for i = 1:Nx
                @inbounds Gᵥ[i,j] += paramKd[i,j,nv,nw] * F[i,j,nw]
            end
        end
    end

    return nothing
end

# 1D
function apply_paramop_supdiag!(G::AbstractArray{<:Number,2},  # output field; G[i,w] is w-component of G at (i,)
                                F::AbstractArray{<:Number,2},  # input field; F[i,w] is w-component of F at (i,)
                                paramKd::AbstractArray{<:Number,3},  # size of last two dimensions is Kf-by-Kf
                                kdiag::Integer;  # 0 ≤ kdiag ≤ Kf-1; index of diagonal of material parameter tensor to set (kdiag = 0: main diagonal)
                                α::Number=1.0,  # scale factor to multiply to result before adding it to G: G += α param F
                                overwrite::Bool=true)  # false to add result to output storage G
    Nx, Kf = size(G)  # Kf: field dimension
    @assert size(paramKd)[end] == size(paramKd)[end-1] == Kf

    overwrite && (G .= 0)
    for nv = 1:Kf  # row index of material parameter tensor
        Gᵥ = @view G[:,:,nv]
        nw = mod1(nv+kdiag, Kf)  # column index of material parameter tensor
        @threads for j = 1:Ny
            for i = 1:Nx
                @inbounds Gᵥ[i,j] += paramKd[i,j,nv,nw] * F[i,j,nw]
            end
        end
    end

    return nothing
end

# function apply_paramop_supdiag!(G::AbstractArray{<:Number,K₊₁},  # output field; in 3D, G[i,j,k,w] is w-component of G at (i,j,k)
#                                 F::AbstractArray{<:Number,K₊₁},  # input field; in 3D, F[i,j,k,w] is w-component of F at (i,j,k)
#                                 paramKd::AbstractArray{<:Number,K₊₂},  # size of last two dimensions is Kf-by-Kf
#                                 kdiag::Integer;  # 0 ≤ kdiag ≤ Kf-1; index of diagonal of material parameter tensor to set (kdiag = 0: main diagonal)
#                                 α::Number=1.0  # scale factor to multiply to result before adding it to G: G += α param F
#                                 overwrite::Bool=true)  # false to add result to output storage G
#     K = K₊₁ - 1
#     @assert K₊₂==K+2
#
#     Kf = size(paramKd, K+1)  # field dimension
#     @assert size(paramKd,K+2)==Kf
#
#     overwrite && (G .= 0)
#     CI = CartesianIndices(size(G)[1:K])
#     for nv = 1:Kf  # row index of material parameter tensor
#         nw = mod1(nv+kdiag, Kf)  # column index of material parameter tensor
#         for ci = CI
#             @inbounds G[ci,nv] += paramKd[ci,nv,nw] * F[ci,nw]
#         end
#     end
#
#     return nothing
# end
