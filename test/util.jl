@testset "util" begin

@testset "sort_ind!" begin
    v = MArray{Tuple{2,2,2}}(rand(8))
    ind = MArray{Tuple{2,2,2}}(collect(1:8))
    @inferred(MaxwellBase.sort_ind!(ind,v))
    @test issorted(v[ind])
    @test @inferred(MaxwellBase.countdiff(ind,v)) == (8,8)
end  # @testset "sort8!"

end  # @testset "util"
