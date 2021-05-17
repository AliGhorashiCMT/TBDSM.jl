using Pkg;
ENV["PYTHON"]=""
Pkg.build("PyCall")
using Test, TBDSM

@testset "Specific Examples" begin
    @test TBDSM.heaviside(1) == 1
    @test TBDSM.heaviside(-1) == 0
end

@testset "Polarization Tests" begin
    K1= [-17.0309799458612, 0.0]
    g1 = impol_2d(1/6, 0, 1,  pb_graphene.monolayer(), spin=4, mesh=30, offset=[0, 0], subsampling=6)
    g2 = graphene_impol(1/6, 0,1, mesh=30, histogram_width=10, subsampling=6)
    @test 100*maximum(replace(x-> x==Inf ? 0 : x, replace(x-> isnan(x) ?  0 : x, (g2-g1)./g1))) < 1e-2 #Less than .01 percent error in polarization calculations
end

@testset "Density of States Tests" begin
    #=
    We check that the density of states method give the correct results for dirac cones. 
    The analytical result is given by 4*A/(4pi^2)*(2pi)*E/v^2. A = 3*sqrt(3)/2*1.42^2
    =#
    es, vs = dos(pb_model(pb_graphene.monolayer(), pb.translational_symmetry()), mesh=300, histogram_width=5)
    esb, vsb = dos(pb_model(pb_graphene.bilayer(), pb.translational_symmetry()), mesh=300, histogram_width=5)
    dosanalytic(x::Real) = 3*sqrt(3)/2*1.42^2*2*pi/(pi^2)*abs(x)/36
    for (e, v) in zip(es, vs)
        ((abs(e)<1) & (v != 0) & (dosanalytic(e) !=0)) && @test abs(100*(dosanalytic(e)-2*v)/(2*v)) < 10 
    end
    #Below we check for up to 5 percent errors in comparing bilayer and monolayer graphene densities of states
    for ϵ in -1.5:0.1:1.5
        println(ϵ)
        ϵ!=0 && (@test 5 > (100*(vsb[argmin(abs.(esb.-ϵ))]-2*vs[argmin(abs.(es.-ϵ))])/(2*vs[argmin(abs.(es.-ϵ))]))) 
    end  
end

@testset "bilayergraphene" begin

    #test that bilayer graphene has correct dispersion as mass term varies
    bilayer = pb_model(bilayer_nointerlayer(0), pb.translational_symmetry())
    graphene = pb_model(pb_graphene.monolayer(), pb.translational_symmetry())
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    bilayer_disp = pb_solver(bilayer).calc_bands(Gamma, K2, M, K1).energy
    monolayer_disp = pb_solver(graphene).calc_bands(Gamma, K2, M, K1).energy
    @test bilayer_disp[:, 1] ≈ monolayer_disp[:, 1]
    @test bilayer_disp[:, 2] ≈ monolayer_disp[:, 1]
    @test bilayer_disp[:, 3] ≈ monolayer_disp[:, 2]
    @test bilayer_disp[:, 4] ≈ monolayer_disp[:, 2]

    #Now make a shifted bilayer dispersion 
    bilayer = pb_model(bilayer_nointerlayer(20), pb.translational_symmetry())
    bilayer_disp = pb_solver(bilayer).calc_bands(Gamma, K2, M, K1).energy

    @test bilayer_disp[:, 1] ≈ monolayer_disp[:, 1]
    @test bilayer_disp[:, 2] ≈ monolayer_disp[:, 2]

    @test bilayer_disp[:, 3] ≈  monolayer_disp[:, 1].+20
    @test bilayer_disp[:, 4] ≈  monolayer_disp[:, 2].+20
end

@testset "ImPolarization Test" begin
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    imps = impol_2d(1/6, 0, 1, bilayer_nointerlayer(.001), spin=4, mesh=200, 
    histogram_width=5, offset=K1, subsampling=3)
    impsmon=impol_2d(1/6, 0, 1, pb_graphene.monolayer(),
    spin=4, mesh=200, histogram_width=5, offset=K1, subsampling=3)
    percentdiff = (2*impsmon-imps)./imps
    @test 100*maximum(abs.(replace(x-> x==Inf ? 0 : x, replace(x-> isnan(x) ?  0 : x, percentdiff)))) < 5 #Less than 5 percent error in polarization calculations
end

@testset "Twisted Bilayer Graphene Model" begin
    b1, b2 = TBDSM.levitov_tbg_model().reciprocal_vectors()
    #Check lengths of reciprocal lattice vectors
    correct_length = 2π/13.4*2/sqrt(3)
    isapprox(sqrt(sum(b1.^2)), correct_length, atol=1e-3)
    isapprox(sqrt(sum(b2.^2)),  correct_length, atol=1e-3)
end

@testset "TMD" begin
    tmd_mo_s2(false)
end