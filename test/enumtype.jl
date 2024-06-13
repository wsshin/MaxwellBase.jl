@testset "enumtype" begin

@testset "instances" begin
    @test instances(AbstractFieldType) == (Efield(), Hfield())
end

@testset "conj" begin
    @test conj(Efield()) == Hfield()
    @test conj(Hfield()) == Efield()
end

@testset "gt4ft" begin
    @test MaxwellBase.gt4ft(Efield(),Efield()) == PrimalGrid()
    @test MaxwellBase.gt4ft(Hfield(),Efield()) == DualGrid()
    @test MaxwellBase.gt4ft(Efield(),Hfield()) == DualGrid()
    @test MaxwellBase.gt4ft(Hfield(),Hfield()) == PrimalGrid()
end

end  # @testset "enumtype"
