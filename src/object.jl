# Below, the symbols of GeometryPrimitives such as KDTree, bounds, surfpt_nearby, normal do
# not have to be re-defined here to be exported.  When GeometryPrimitives is imported to
# MaxwellBase by "using GeometryPrimitives", the symbols of GeometryPrimitives are already
# inside the namespace of MaxwellBase, so we simply need to export them to use outside
# MaxwellBase without qualification.
#
# One caveat, though.  In this file, bounds is defined for Interval, so if we simply export
# bounds, then only this is exported.  So, when defining bounds for Interval, we have to
# "extend" bounds for Shape by defining it as GeometryPrimitives.bounds(::Interval) = ....
# Then, exporting bounds exports this whole collection of bounds, both for Shape and Interval.
export Object
export shape, matparam, max∆l, sub_matprm4pind, add_obj!, periodize  #, surfpt_nearby, normal
# export lsf, bound_, L_, center_, dist2bound, bound_contains, ∆lmax, sphere, transform,
#     surfnormal, surfpoint  # functions
# import Base:size, getindex, contains, isless, union, intersect

mutable struct Object{K,S<:Shape{K}}  # K: dimension of space
    shp::S
    mat::Material
    ∆lmax::SVector{K,Float64}
    Object{K,S}(shp, mat, ∆lmax) where {K,S} = new(shp, mat, ∆lmax)
end

Object(shp::S, mat::Material, ∆lmax::SVector{K}=fill(Inf,SVector{K})) where {K,S<:Shape{K}} =
    Object{K,S}(shp, mat, ∆lmax)
Object(shp::Shape{K}, mat::Material, ∆lmax::AbstractVector) where {K} = Object(shp, mat, SVector{K}(∆lmax))
Object(shp::Shape{K}, mat::Material, ∆lmax::Real) where {K} = Object(shp, mat, SVector(ntuple(k->∆lmax, Val(K))))

# Create the user interface such that maxwellsys.add(shape, material, ∆lmax) uses this function.
# setmat!(obj::Object, m::Material) = (obj.mat = m; obj)
# setmax∆l!(obj::Object{K}, ∆lmax::AbstractVector) where {K} = (obj.∆lmax = SVector{K}(∆lmax); obj)
# setmax∆l!(obj::Object{K}, ∆lmax::Number) where {K} = setmax∆l!(obj, SVector(ntuple(k->∆lmax, Val(K))))

shape(obj::Object) = obj.shp
matparam(obj::Object, ft::AbstractFieldType) = matparam(obj.mat, ft)
max∆l(obj::Object) = obj.∆lmax

# Relationships between matprm4pind, pind4oind, shp4oind
#
# - matprm4pind is a vector of material parameter tensors.  matprm4pind[n] is a 3×3 material
# parameter tensor for the parameter index n.
#
# - pind4oind is a map from the object index to the parameter index.
#
# - shp4oind is a map from the object index to the shape.

# Returns a vector of submatrices of material parameter tensors.
# We return a vector of submatrices rather than a vector of submatrix views, because
# assign_param!() requires the entries of its argument matprm4pind to be SMatrix rather than
# views of SMatrix.
sub_matprm4pind(
    matprm4pind::AbstractVector{<:SMatrix{Kf,Kf,<:Number,Kf²}},
    inds::SVector{<:Any,<:Integer}
) where {Kf,Kf²} =
    [mp[inds,inds] for mp = matprm4pind]  # mp[inds,inds]: SMatrix{K,K}

# Consider using resize! on shp4oind.
function add_obj!(
    shp4oind::AbstractVector{Shape{K,K²}},
    pind4oind::NTuple{2,AbstractVector{ParamInd}},
    matprm4pind::NTuple{2,AbstractVector{<:SMatrix{3,3,<:Number,9}}},
    objs::AbstractVector{Object{K,<:Shape{K,K²}}}
) where {K,K²}
    for obj = objs
        add_obj!(shp4oind, pind4oind, matprm4pind, obj)
    end
end

# Consider using resize! on shp4oind.
function add_obj!(
    shp4oind::AbstractVector{Shape{K,K²}},
    pind4oind::NTuple{2,AbstractVector{ParamInd}},
    matprm4pind::NTuple{2,AbstractVector{<:SMatrix{3,3,<:Number,9}}},
    objs::Object{K,<:Shape{K,K²}}...
) where {K,K²}
    for obj = objs
        add_obj!(shp4oind, pind4oind, matprm4pind, obj)
    end
end

# Populate `shp4oind`, `pind4oind`, `matprm4pind` with a given object.  See above for the
# meaning of `shp4oind`, `pind4oind`, `matprm4pind`.  In doing so, figure out the unique
# indices for the material parameters and use them as pind in `pind4oind` and `matprm4pind`.
#
# When I put an periodic array of objects, some object could be across a boundary.  In that
# case, the rule is to put two different objects at the negative-end and positive-end
# boundaries.  Make sure, though, the two objects have the same surface normal at the
# boundary points, so that the surface normal calculation is not ambiguous there.
#
# Note that matprm4pind stores the same material parameter "value" only once.  Therefore,
# if two objects are created with different Material instances with the same material
# parameter tensor, then the two objects' materials are assigned with the same pind.
function add_obj!(
    shp4oind::AbstractVector{Shape{K,K²}},  # initially empty vector
    pind4oind::NTuple{2,AbstractVector{ParamInd}},  # tuple of initially empty vectors
    matprm4pind::NTuple{2,AbstractVector{<:SMatrix{3,3,<:Number,9}}},  # tuple of initially empty vectors
    obj::Object{K,<:Shape{K,K²}}
) where {K,K²}
    push!(shp4oind, shape(obj))  # append obj (for potential use of obj4oind with KDTree, must use pushfirst! to prepend)

    # Assign the material parameter indices to obj.
    # If the material is aleady used in some object that was previously added, then use the
    # index already assigned to that material.
    for ft = instances(AbstractFieldType)
        mp, on2pn, pn2mp = matparam(obj,ft), pind4oind[ft], matprm4pind[ft]
        pind = findlast(x -> x==mp, pn2mp)  # number if mp is already in pn2mp
        pind===nothing && (push!(pn2mp, mp); pind = length(pn2mp))  # update pind if mp is not in pn2mp
        push!(on2pn, pind)
    end

    return nothing
end

function GeometryPrimitives.periodize(obj::Object{K}, A::AbstractMatrix, ∆range::Shape{K}) where {K}
    shp_array = periodize(obj.shp, A, ∆range)
    N = length(shp_array)
    obj_array = Vector{Object{K}}(N)
    for n = 1:N
        obj_array[n] = Object(shp_array[n], obj.mat, obj.∆lmax)
    end

    return obj_array
end


# The following includes operation on shapes.  These will need to be defined in GeometryPrimitives.
# Then, I will need to define the same functions for Object that delegate the operations to
# Shape's functions.
#
# See also test/object.jl


# # distance to the boundary.  If I want to get 0 for internal points, compare
# # dist(l,c) = (distance between l and center), and if it is less than L_(i) / 2,
# # return zero; otherwise return dist(l,c).
# @inline dist2bound(i::Interval1D, l::Real) = (a1 = abs(i.bound[LowerBound()]-l)) < (a2 = abs(i.bound[UpperBound()]-l)) ? a1 : a2

# # Functions common to all Shape's
#

# # contains
# # Sometimes, we only have level-set functions rather than Shape instances.
# @inline contains(lsf::Function, l::NTuple{3,Real}, isinclusive::Bool=true) = isinclusive ? lsf(l) ≥ 0 : lsf(l) > 0
# @inline contains(s::Shape, l::NTuple{3,Real}, isinclusive::Bool=true) = contains(lsf(s), l, isinclusive)
#
# # Other functions
# @inline bound_(s::Shape) = bound_(s.cbox)
# @inline bound_(s::Shape, w::Axis) = bound_(s.cbox, w)
# @inline bound_(s::Shape, w::Axis, bt::AbstractBoundType) = bound_(s.cbox, w, bt)
# @inline bound_contains(s::Shape, l::NTuple{3,Real}, isinclusive::Bool=true) = contains(s.cbox, l, isinclusive)
# @inline L_(s::Shape) = L_(s.cbox)  # scalar
# @inline center_(s::Shape) = center_(s.cbox)  # scalar
#
# # Level-set function manipulation
# @inline transform(lsf::Function, σ::NTuple{3,Real}, ∆::NTuple{3,Real}) = l -> lsf(σ .* (l .- ∆))  # σ[w] = 1 for shift, σ[w] = -1 for flip
#
# @inline union(lsf::Function...) = l -> maximum(map(f->f(l), lsf))
# @inline union(lsf1::Function, lsf2::Function) = lsf1==lsf2 ? lsf1 : l->max(lsf1(l), lsf2(l))
#
# @inline intersect(lsf::Function...) = l -> minimum(map(f->f(l), lsf))
#
# # Numerical gradient, and closest boundary point to a given point (needed for subpixel
# # smoothing of material parameters)
#
# # Calculates the direction normal to the curved surface crossing a voxel by numerical gradient
# # of the level-set function describing the surface.  The gradient is evaluated at a given point
# # inside the voxel.
# function surfnormal(lsf::Function, r::NTuple{3,Float64}, vxl::NTuple{3,NTuple{2,Float64}})
#     # println("lsf((0,0,0)) = $(lsf((0,0,0))),
#     #     lsf((-1,0,0)) = $(lsf((-1,0,0))), lsf((1,0,0)) = $(lsf((1,0,0))),
#     #     lsf((0,-1,0)) = $(lsf((0,-1,0))), lsf((0,1,0)) = $(lsf((0,1,0))),
#     #     lsf((0,0,-1)) = $(lsf((0,0,-1))), lsf((0,0,1)) = $(lsf((0,0,1))),
#     #     lsf((-1,-1,-1)) = $(lsf((-1,-1,-1))), lsf((1,-1,-1)) = $(lsf((1,-1,-1))),
#     #     lsf((-1,1,-1)) = $(lsf((-1,1,-1))), lsf((1,1,-1)) = $(lsf((1,1,-1))),
#     #     lsf((-1,-1,1)) = $(lsf((-1,-1,1))), lsf((1,-1,1)) = $(lsf((1,-1,1))),
#     #     lsf((-1,1,1)) = $(lsf((-1,1,1))), lsf((1,1,1)) = $(lsf((1,1,1))),
#     #     ")
#     lsf2(s) = 1-hypot(s[1], s[2], s[3])
#     # println("lsf($r) = $(lsf(r))")
#     # println("lsf2($r) = $(lsf2(r))")
#     # println("r = $r, vxl = $vxl")
#     α = 5e-4  # 1e-3 / 2
#
#     ∆x, ∆y, ∆z = α * (vxl[nX][UpperBound()] - vxl[nX][LowerBound()]), α * (vxl[nY][UpperBound()] - vxl[nY][LowerBound()]), α * (vxl[nZ][UpperBound()] - vxl[nZ][LowerBound()])
#     # ∆x, ∆y, ∆z = α .* (t_ind(vxl, UpperBound(), UpperBound(), UpperBound()) .- t_ind(vxl, LowerBound(), LowerBound(), LowerBound()))
#     rx, ry, rz = r
#
#     # println("(∆x, ∆y, ∆z) = $((∆x, ∆y, ∆z))")
#     # println("(rx, ry, rz) = $((rx, ry, rz))")
#     # # println("lsf(($rx, $ry, $rz)) = $(lsf((rx,ry,rz)))")
#     fxn::Float64, fxp::Float64 = lsf((rx-∆x, ry, rz)), lsf((rx+∆x, ry, rz))
#     fyn::Float64, fyp::Float64 = lsf((rx, ry-∆y, rz)), lsf((rx, ry+∆y, rz))
#     fzn::Float64, fzp::Float64 = lsf((rx, ry, rz-∆z)), lsf((rx, ry, rz+∆z))
#     # println("((fxn, fxp), (fyn, fyp), (fzn, fzp)) = $(((fxn, fxp), (fyn, fyp), (fzn, fzp)))")
#
#     fxn2::Float64, fxp2::Float64 = lsf2((rx-∆x, ry, rz)), lsf2((rx+∆x, ry, rz))
#     fyn2::Float64, fyp2::Float64 = lsf2((rx, ry-∆y, rz)), lsf2((rx, ry+∆y, rz))
#     fzn2::Float64, fzp2::Float64 = lsf2((rx, ry, rz-∆z)), lsf2((rx, ry, rz+∆z))
#     # println("((fxn2, fxp2), (fyn2, fyp2), (fzn2, fzp2)) = $(((fxn2, fxp2), (fyn2, fyp2), (fzn2, fzp2)))")
#
#     neg∂f∂x, neg∂f∂y, neg∂f∂z = (fxn - fxp) / ∆x, (fyn - fyp) / ∆y, (fzn - fzp) / ∆z
#
#     return neg∂f∂x, neg∂f∂y, neg∂f∂z  # outward normal is negative gradient
# end
#
# # Find the point on a curved surface by moving from a point along a direction.
# function surfpoint(lsf::Function, r::NTuple{3,Float64}, v::NTuple{3,Float64})
#     # println("r = $r, v = $v")
#     f(t) = lsf(r .+ t.*v)
#     tsol, isconverged = newtsol(0., f)
#
#     p = (r[nX] + tsol*v[nX], r[nY] + tsol*v[nY], r[nZ] + tsol*v[nZ])
#
#     return p, isconverged
#     # return r .+ tsol.*v
# end
