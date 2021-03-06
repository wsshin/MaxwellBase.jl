@testset "material" begin

@testset "tensorize" begin
    @test MaxwellBase.tensorize(3, Val(2)) === SMatrix{2,2,ComplexF64,4}(3,0,0,3)
    @test MaxwellBase.tensorize([2,4], Val(2)) === SMatrix{2,2,ComplexF64,4}(2,0,0,4)
    @test MaxwellBase.tensorize([1 3; 2 4], Val(2)) === SMatrix{2,2,ComplexF64,4}(1,2,3,4)
end

@testset "Material" begin
    m_nn = Material{3,3}("(num,num)", ε = 2, μ = 3)
    m_vv = Material{3,3}("(vec,vec)", ε = [2.1,2.2,2.3], μ = [3.1,3.2,3.3])
    m_mm = Material{3,3}("(mat,mat)", ε = [2 2 2; 2 2 2; 2 2 2], μ = [3 3 3; 3 3 3; 3 3 3])

    @test matparam(m_nn,EE) == diagm(0=>[2,2,2]) && matparam(m_nn,HH) == diagm(0=>[3,3,3]) && string(m_nn) == "(num,num)"
    @test matparam(m_vv,EE) == diagm(0=>[2.1,2.2,2.3]) && matparam(m_vv,HH) == diagm(0=>[3.1,3.2,3.3]) && string(m_vv) == "(vec,vec)"
    @test matparam(m_mm,EE) == fill(2, (3,3)) && matparam(m_mm,HH) == fill(3, (3,3)) && string(m_mm) == "(mat,mat)"
end

@testset "τ_trans" begin
    # Check if τ_trans for a 2×2 matrix is τ_trans for the 3×3 version of the matrix.  The
    # 3×3 version needs to have the original 2×2 matrix as the top-left 2×2 block, and the
    # third row and column has the (3,3) entry as the only nonzero entry.
    ε2 = @SMatrix rand(ComplexF64, 2, 2)
    ε3temp = zeros(ComplexF64, 3, 3)
    ε3temp[1:2,1:2] .= ε2
    ε3temp[3,3] = rand(ComplexF64)
    ε3 = MaxwellBase.S²ComplexF3(ε3temp)

    τ2 = MaxwellBase.τ_trans(ε2)
    τ3 = MaxwellBase.τ_trans(ε3)
    @test τ3[1:2,1:2] == τ2
end

end  # @testset "material"
