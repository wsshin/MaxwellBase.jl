@testset "object" begin

@testset "Object" begin
    mat = Material("material", ε=rand(), μ=rand())

    @testset "Cuboid" begin
        ba = sort(rand(3,2), dims=2)  # array
        b = (ba[:,1], ba[:,2])

        box_mat = Object(Cuboid(b), mat)  # new convenience constructor

        @test shape(box_mat) == Cuboid(b)
        @test max∆l(box_mat) == fill(Inf,3)
        @test matparam(box_mat,EE)==mat.param[nE] && matparam(box_mat,HH)==mat.param[nH]
    end  # @testset "Cuboid"
end  # @testset "Object"

# @testset "periodize" begin
#     # Square lattice
#     Si = Material("Si", ε = 12)
#     gtₑ = PrimalGrid()
#     eSi = EncodedMaterial(gtₑ, Si)
#     c_eSi = Object(Cylinder([0,0,0], 1, 5, [0,0,1]), eSi)
#
#     b = Cuboid([0,0,0], [10,10,5])
#     A = [1 0 0; 0 1 0; 0 0 5]'
#     obj_array = periodize(c_eSi, A, b)
#
#     @test length(obj_array) == 11^2
#     @test all(map(x->(x==max∆l(c_eSi)), max∆l.(obj_array)))
#     @test all(map(x->(x==matparam(c_eSi,PrimalGrid())), matparam.(obj_array,PrimalGrid())))
#     @test all(map(x->(x==matparam(c_eSi,DualGrid())), matparam.(obj_array,DualGrid())))
#
#     obj4oind = Object{3}[]
#     matprm4pind = (S²ComplexF3[], S²ComplexF3[])
#     add_obj!(obj4oind, matprm4pind, obj_array)
# end  # @testset "periodize"


# @testset "union" begin
#     b1 = Cuboid(((0,1),(0,1),(0,1)))
#     b2 = Cuboid(((-1,0),(-1,0),(-1,0)))
#     s = sphere((0,0,0), 1)
#     lsf_union = union(lsf(b1), lsf(s), lsf(b2))
#     @test all(contains.(lsf_union, [(1,1,1), (-1,-1,-1), (-1/√3,-1/√3,1/√3), (1/√3,1/√3,-1/√3)]))
# end
#
# @testset "intersect" begin
#     b1 = Cuboid(((-1,2),(-1,2),(-1,2)))
#     b2 = Cuboid(((-2,1),(-2,1),(-2,1)))
#     lsf_intersect = intersect(lsf(b1), lsf(b2))
#     @test all(contains.(lsf_intersect, [(x,y,z) for x = -1:1, y = -1:1, z = -1:1]))
#     @test all(!contains.(lsf_intersect, [(x,y,z) for x = (-2,2), y = (-2,2), z = (-2,2)]))
# end
#
# @testset "flip" begin
#     b = Cuboid(((0,1),(0,1),(0,1)))
#     lsf0 = lsf(b)
#     lsf_x = flip(lsf0, X̂, 0)
#     lsf_y = flip(lsf0, Ŷ, 0)
#     lsf_z = flip(lsf0, Ẑ, 0)
#
#     @test all(contains.(lsf0, [(x,y,z) for x = 0:1, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_x, [(x,y,z) for x = -1:0, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_y, [(x,y,z) for x = 0:1, y = -1:0, z = 0:1]))
#     @test all(contains.(lsf_z, [(x,y,z) for x = 0:1, y = 0:1, z = -1:0]))
# end
#
# @testset "shift" begin
#     b = Cuboid(((0,1),(0,1),(0,1)))
#     lsf0 = lsf(b)
#     lsf_x = shift(lsf0, X̂, -1)
#     lsf_y = shift(lsf0, Ŷ, -1)
#     lsf_z = shift(lsf0, Ẑ, -1)
#
#     @test all(contains.(lsf0, [(x,y,z) for x = 0:1, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_x, [(x,y,z) for x = -1:0, y = 0:1, z = 0:1]))
#     @test all(contains.(lsf_y, [(x,y,z) for x = 0:1, y = -1:0, z = 0:1]))
#     @test all(contains.(lsf_z, [(x,y,z) for x = 0:1, y = 0:1, z = -1:0]))
# end

end  # @testset "object"
