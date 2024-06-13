@testset "param" begin

# Need to test symmetry of param3dmat for different boundary conditions.
@testset "symmetry of create_paramop, uniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Cuboid(g3.bounds), vac)
    obj_diel = Object(Cuboid([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Ball([0,0,0], 1), diel)

    # Add objects.
    shp4oind = Shape3[]
    εind4oind = ParamInd[]
    μind4oind = ParamInd[]
    ε4εind = SMatrix{3,3,ComplexF64,9}[]
    μ4μind = SMatrix{3,3,ComplexF64,9}[]

    add_obj!(shp4oind, (εind4oind,μind4oind), (ε4εind,μ4μind), dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    N = g3.N
    ε3d = create_param_array(ComplexF64, N)
    oind3d_εxx = create_oind_array(N)
    oind3d_εyy = create_oind_array(N)
    oind3d_εzz = create_oind_array(N)
    oind3d_εoo = create_oind_array(N)

    μ3d = create_param_array(ComplexF64, N)
    oind3d_μxx = create_oind_array(N)
    oind3d_μyy = create_oind_array(N)
    oind3d_μzz = create_oind_array(N)
    oind3d_μoo = create_oind_array(N)

    boundft = SVector(Efield(),Efield(),Efield())
    assign_param!(ε3d, (oind3d_μxx,oind3d_μyy,oind3d_μzz), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(ε3d, tuple(oind3d_μoo), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(μ3d, (oind3d_εxx,oind3d_εyy,oind3d_εzz), shp4oind, μind4oind, μ4μind, ft2gt.(Hfield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(μ3d, tuple(oind3d_εoo), shp4oind, μind4oind, μ4μind, ft2gt.(Hfield(),boundft), g3.ghosted.τlg, g3.isbloch)

    # Perform smoothing.
    smooth_param!(ε3d, (oind3d_εxx,oind3d_εyy,oind3d_εzz), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.l, g3.ghosted.lg, g3.σ, g3.ghosted.∆τ)
    smooth_param!(ε3d, tuple(oind3d_εoo), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.l, g3.ghosted.lg, g3.σ, g3.ghosted.∆τ)

    ∆l = g3.∆l[DualGrid()]
    ∆l′ = g3.∆l[PrimalGrid()]
    ∆l′⁻¹ = map(x->inv.(x), ∆l′)

    Mε = create_paramop(ε3d, boundft.!=Efield(), ∆l, ∆l′⁻¹, g3.isbloch, cmpfirst=true)

    # @test issymmetric(Mε)
    @test Mε ≈ transpose(Mε)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mε3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "create_paramop"

@testset "symmetry of create_paramop, nonuniform grid" begin
    # Create a grid.
    isbloch = [true, true, true]
    # isbloch = [true, false, false]
    lprim = ([-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5], [-1.5, -0.4, 0.3, 1.5])
    g3 = Grid(lprim, isbloch)
    N = g3.N

    # Create materials.
    εvac = 1.0
    vac = Material("Vacuum", ε=εvac)

    εdiel = 2.0
    diel = Material("Dielectric", ε=εdiel)

    # Create objects.
    dom_vac = Object(Cuboid(g3.bounds), vac)
    obj_diel = Object(Cuboid([0,0,0], [1,1,1]), diel)
    # obj_diel = Object(Ball([0,0,0], 1), diel)

    # Add objects.
    shp4oind = Shape3[]
    εind4oind = ParamInd[]
    μind4oind = ParamInd[]
    ε4εind = SMatrix{3,3,ComplexF64,9}[]
    μ4μind = SMatrix{3,3,ComplexF64,9}[]

    add_obj!(shp4oind, (εind4oind,μind4oind), (ε4εind,μ4μind), dom_vac, obj_diel)

    # Construct arguments and call assign_param!.
    N = g3.N
    ε3d = create_param_array(ComplexF64, N)
    oind3d_εxx = create_oind_array(N)
    oind3d_εyy = create_oind_array(N)
    oind3d_εzz = create_oind_array(N)
    oind3d_εoo = create_oind_array(N)

    μ3d = create_param_array(ComplexF64, N)
    oind3d_μxx = create_oind_array(N)
    oind3d_μyy = create_oind_array(N)
    oind3d_μzz = create_oind_array(N)
    oind3d_μoo = create_oind_array(N)

    boundft = SVector(Efield(),Efield(),Efield())
    assign_param!(ε3d, (oind3d_μxx,oind3d_μyy,oind3d_μzz), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(ε3d, tuple(oind3d_μoo), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(μ3d, (oind3d_εxx,oind3d_εyy,oind3d_εzz), shp4oind, μind4oind, μ4μind, ft2gt.(Hfield(),boundft), g3.ghosted.τlg, g3.isbloch)
    assign_param!(μ3d, tuple(oind3d_εoo), shp4oind, μind4oind, μ4μind, ft2gt.(Hfield(),boundft), g3.ghosted.τlg, g3.isbloch)

    # Perform smoothing.
    smooth_param!(ε3d, (oind3d_εxx,oind3d_εyy,oind3d_εzz), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.l, g3.ghosted.lg, g3.σ, g3.ghosted.∆τ)
    smooth_param!(ε3d, tuple(oind3d_εoo), shp4oind, εind4oind, ε4εind, ft2gt.(Efield(),boundft), g3.l, g3.ghosted.lg, g3.σ, g3.ghosted.∆τ)

    ∆l = g3.∆l[DualGrid()]
    ∆l′ = g3.∆l[PrimalGrid()]
    ∆l′⁻¹ = map(x->inv.(x), ∆l′)

    Mε = create_paramop(ε3d, boundft.!=Efield(), ∆l, ∆l′⁻¹, g3.isbloch, cmpfirst=true)

    # Ml
    ∆lprim = g3.∆l[PrimalGrid()]
    ∆ldual = g3.∆l[DualGrid()]

    Nx, Ny, Nz = N
    nX, nY, nZ = 1, 2, 3
    ∆X = repeat(∆ldual[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(reshape(∆ldual[nY], (1,Ny,1)), outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆ldual[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆L = [reshape(∆X,1,:); reshape(∆Y,1,:); reshape(∆Z,1,:)]
    ∆l = ∆L[:]
    Ml⁻¹ = sparse(Diagonal(1 ./ ∆l))

    # Ma
    ∆X = repeat(∆lprim[nX], outer=(1,Ny,Nz))
    ∆Y = repeat(reshape(∆lprim[nY], (1,Ny,1)), outer=(Nx,1,Nz))
    ∆Z = repeat(reshape(∆lprim[nZ], (1,1,Nz)), outer=(Nx,Ny,1))

    ∆YZ = ∆Y .* ∆Z
    ∆ZX = ∆Z .* ∆X
    ∆XY = ∆X .* ∆Y

    ∆A = [reshape(∆YZ,1,:); reshape(∆ZX,1,:); reshape(∆XY,1,:)]
    ∆a = ∆A[:]
    Ma = sparse(Diagonal(∆a))

    Mε_sym = Ma * Mε * Ml⁻¹

    @test Mε_sym ≈ transpose(Mε_sym)

    # Code I used to reveal the matrix structure while debugging:
    # (i,j) = (nX,nY); full(real.(Mε3))[27(i-1)+1:27i, 27(j-1)+1:27j]
end  # @testset "create_paramop"

# Need to test if PML still works after the field-averaging operator applied to the material
# parameters, because effectively it makes the material inhomogeneous along the direction
# normal to the PML interface.

end  # @testset "param"
