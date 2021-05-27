@testset "phys" begin

L₀ = 1e-9
λ = 1550  # 1550L₀
osc = Oscillation(1550, L₀)

@test in_L₀(osc) ≈ λ
@test in_ω₀(osc) ≈ 2π / λ
@test in_eV(osc) ≈ (ℎ*c₀) / (λ*L₀)

end  # @testset "phys"
