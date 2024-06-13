@testset "enumtype" begin

@testset "instances" begin
    @test instances(AbstractFieldType) == (Efield(), Hfield())
end

@testset "conj" begin
    @test conj(Efield()) == Hfield()
    @test conj(Hfield()) == Efield()
end

@testset "ft2gt" begin
    @test MaxwellBase.ft2gt(Efield(),Efield()) == PrimalGrid()
    @test MaxwellBase.ft2gt(Hfield(),Efield()) == DualGrid()
    @test MaxwellBase.ft2gt(Efield(),Hfield()) == DualGrid()
    @test MaxwellBase.ft2gt(Hfield(),Hfield()) == PrimalGrid()
end

end  # @testset "enumtype"
