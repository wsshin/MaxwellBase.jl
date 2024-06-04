export Material
export matparam, kottke_avg_param

# Note that the output of tensorize is always a matrix: for K = 1, it is a 1×1 matrix.
tensorize(x::Number) = tensorize(fill(x,SVector{3}))  # SComplexF{3}(x,x,x): vector
tensorize(v::AbstractVector{<:Number}) = diagm(Val(0)=>SVector{3}(v))
tensorize(m::AbstractMatrix{<:Number}) = SMatrix{3,3}(m)

struct Material
    name::String
    param::NTuple{2,SMatrix{3,3,<:Number,9}}  # (material parameter interacting with E, and H)
end

Material(name::String, param::NTuple{2,AbstractMatParam}) = Material(name, tensorize.(param))
Material(name::String; ε::AbstractMatParam=1, μ::AbstractMatParam=1) = Material(name, (ε,μ))  # ε, μ: keyword arguments
# Consider implementing Material2 and Material3 for 2D and 3D solvers.

Base.string(m::Material) = m.name
matparam(m::Material, ft::FieldType) = m.param[Int(ft)]

kottke_avg_param(
    param1::AbstractMatrix{<:Number},
    param2::AbstractMatrix{<:Number},
    n12::AbstractVector{<:Real},
    rvol1::Real
) = (
    K = length(n12);
    kottke_avg_param(SMatrix{K,K}(param1), SMatrix{K,K}(param2), SVector{K}(n12), rvol1)
)

# Implement the averaging scheme between two local material parameter tensors developed in
# the paper by Kottke, Farjadpour, Johnson entitled "Perturbation theory for anisotropic
# dielectric interfaces, and application to subpixel smoothing of discretized numerical
# methods", Physical Review E 77 (2008): 036611.
function kottke_avg_param(
    param1::SMatrix{K,K,<:Number,K²},
    param2::SMatrix{K,K,<:Number,K²},
    n12::SVector{K,Float64},
    rvol1::Real
) where {K,K²}
    @argcheck K² == K^2
    Scomp = rand(SMatrix{K,K-1})  # directions complementary to n12; works even for K = 1
    Stemp = [n12 Scomp]  # SMatrix{K,K}
    S = qr(Stemp).Q  # nonallocating; 1st column is normalized n12

    τ1 = τ_trans(transpose(S) * param1 * S)  # express param1 in S coordinates, and apply τ transform
    τ2 = τ_trans(transpose(S) * param2 * S)  # express param2 in S coordinates, and apply τ transform

    τavg = τ1 .* rvol1 + τ2 .* (1-rvol1)  # volume-weighted average

    return S * τ⁻¹_trans(τavg) * transpose(S)  # apply τ⁻¹ and transform back to global coordinates
end

# Kottke's averaging scheme reduced for the field dimension orthogonal to the shape
# dimension (such as Hz in the TE mode).  In that case, the field is always parallel to the
# shape boundaries, so the averaging scheme reduces to simple arithmetic averaging.  For
# example for the case with Shape{2} and Kf = 1 (like Hz in the TE mode), the material
# parameter is an 1×1 tensor (= scalar) and therefore the situation is the same as averaging
# isotropic ε between E-field voxels, which is arithmetic averaging.  The only other case is
# the case with Shape{1} and Kf = 1  (like Ex in the slapes in the z-direction).  We can
# easily prove that in these cases again the averaging reduces to simple arithmetic
# averaging between material parameter tensors.  See Agenda > MaxwellBase > Feb/22/2020.
kottke_avg_param(param1::AbstractMatrix{<:Number}, param2::AbstractMatrix{<:Number}, rvol1::Real) =
    param1 .* rvol1 + param2 .* (1-rvol1)

# Equation (4) of the paper by Kottke et al.  Need to verify for K = 1 and 2.
function τ_trans(ε::SMatrix{3,3,<:Number})
    ε₁₁, ε₂₁, ε₃₁, ε₁₂, ε₂₂, ε₃₂, ε₁₃, ε₂₃, ε₃₃ = ε

    return SMatrix{3,3}(
        -1/ε₁₁, ε₂₁/ε₁₁, ε₃₁/ε₁₁,
        ε₁₂/ε₁₁, ε₂₂ - ε₂₁*ε₁₂/ε₁₁, ε₃₂ - ε₃₁*ε₁₂/ε₁₁,
        ε₁₃/ε₁₁, ε₂₃ - ε₂₁*ε₁₃/ε₁₁, ε₃₃ - ε₃₁*ε₁₃/ε₁₁
    )
    # t = zeros(3,3)
    # t[2:3, 2:3] = s[2:3, 2:3]
    # s11 = s[1,1]
    # s[1,1] = -1
    #
    # t = t - (s(:,1) * s(1,:)) ./ s11
    # return t
end

function τ_trans(ε::SMatrix{2,2,<:Number})
    ε₁₁, ε₂₁, ε₁₂, ε₂₂ = ε

    return SMatrix{2,2}(
        -1/ε₁₁, ε₂₁/ε₁₁,
        ε₁₂/ε₁₁, ε₂₂ - ε₂₁*ε₁₂/ε₁₁
    )
end

τ_trans(ε::SMatrix{1,1,<:Number}) = -1 ./ ε  # SMatrix{1,1,<:Number}

# function τ_trans(ε::SMatrix{K,K,ComplexF64}) where {K}
#     τ_temp1 = @SMatrix [ε[i,j] for i = 2:K, j = 2:K]  # cannot use type parameter K here: error is generated
#     τ_temp2 = [zeros(SMatrix{1,K-1,ComplexF64}); τ_temp1]
#     τ_temp = [zeros(SVector{K,ComplexF64}) τ_temp2]
#
#     σcol_temp = @SVector [ε[i,1] for i = 2:K]
#     σcol = [SVector{1,ComplexF64}(-1); σcol_temp]
#
#     σrow_temp = @SVector [ε[1,j] for j = 2:K]
#     σrow = [SVector{1,ComplexF64}(-1); σrow_temp]
#
#     τ = τ_temp - (σcol * transpose(σrow)) ./ ε[1,1]
#
#     return τ
# end

# Equation (23) of the paper by Kottke et al.
function τ⁻¹_trans(τ::SMatrix{3,3,<:Number})
    τ₁₁, τ₂₁, τ₃₁, τ₁₂, τ₂₂, τ₃₂, τ₁₃, τ₂₃, τ₃₃ = τ

    return SMatrix{3,3}(
        -1/τ₁₁, -τ₂₁/τ₁₁, -τ₃₁/τ₁₁,
        -τ₁₂/τ₁₁, τ₂₂ - τ₂₁*τ₁₂/τ₁₁, τ₃₂ - τ₃₁*τ₁₂/τ₁₁,
        -τ₁₃/τ₁₁, τ₂₃ - τ₂₁*τ₁₃/τ₁₁, τ₃₃ - τ₃₁*τ₁₃/τ₁₁
    )

    # s = zeros(3,3)
    # s[2:3, 2:3] = t[2:3, 2:3]
    # t11 = t[1,1]
    # t[1,1] = 1
    #
    # s = s - (t(:,1) * t(1,:)) ./ t11
    # return s
end

function τ⁻¹_trans(τ::SMatrix{2,2,<:Number})
    τ₁₁, τ₂₁, τ₁₂, τ₂₂ = τ

    return SMatrix{2,2}(
        -1/τ₁₁, -τ₂₁/τ₁₁,
        -τ₁₂/τ₁₁, τ₂₂ - τ₂₁*τ₁₂/τ₁₁
    )
end

τ⁻¹_trans(τ::SMatrix{1,1,<:Number}) = -1 ./ τ  # SMatrix{1,1,<:Number}
