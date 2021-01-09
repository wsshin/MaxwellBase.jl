@testset "field" begin

Nin = (5,6,7)
Min = prod(Nin)
ψ = randn(Nin)
w = instances(Axis)[rand(1:3)]

@testset "primal ghost appender" begin
    Nout = Nin .+ (instances(Axis).==w)
    f = randn()  # factor to multiply

    # Check generation of ghost fields.
    ψghost = create_ghosts_prim(ψ, w, f)
    nw = Int(w)
    outrange = (v->(Ninv = Nin[Int(v)]; v≠w ? (1:Ninv) : (1:1))).(XYZ)
    @test ψghost == ψ[outrange...] .* f


    # Ψ = append_ghosts_prim(ψ, w, f)
    # Aapp = create_ghostappender_prim(Nin, w, f)
    # Arem = create_ghostremover_prim(Nin, w)
    # @test Ψ ≈ reshape(Aapp * ψ[:], Nout)
    # @test Matrix{Float64}(I, Min, Min) ≈ Arem * Aapp
end

@testset "dual ghost appender" begin
    Nout = Nin .+ 2.*(instances(Axis).==w)
    f1, fN = randn(), randn()  # factors to multiply
    isbloch = rand(Bool)

    # Check generation of ghost fields.
    ψghost1, ψghostN = create_ghosts_dual(ψ, w, isbloch, (f1, fN))
    nw = Int(w)
    Ninw = Nin[nw]
    outrange1 = (v->(Ninv = Nin[Int(v)]; v≠w ? (1:Ninv) : (1:1))).(XYZ)
    outrangeN = (v->(Ninv = Nin[Int(v)]; v≠w ? (1:Ninv) : (Ninw:Ninw))).(XYZ)
    if isbloch
        @test ψghost1 == ψ[outrangeN...] .* fN
        @test ψghostN == ψ[outrange1...] .* f1
    else
        @test ψghost1 == ψ[outrange1...] .* f1
        @test ψghostN == ψ[outrangeN...] .* fN
    end

    # Ψ = append_ghosts_dual(ψ, w, true, (f₁, f₂))
    # Aapp = create_ghostappender_dual(Nin, w, true, (f₁, f₂))
    # @test Ψ ≈ reshape(Aapp * ψ[:], Nout)
end

end  # @testset "field"
