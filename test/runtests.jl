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
    g1 = impol_2d(1/6, 0, 1,  pb_graphene.monolayer(), spin=4, mesh=200, offset=K1, subsampling=3)
    g2 = graphene_impol(1/6, 0,1, mesh=200, histogram_width=10)
    @test 100*maximum(replace(x-> x==Inf ? 0 : x, replace(x-> isnan(x) ?  0 : x, (g2-g1)./g1))) < 1e-2 #Less than .01 percent error in polarization calculations
end

@testset "Density of States Tests" begin
    #=
    We check that the density of states method give the correct results for dirac cones. 
    The analytical result is given by 4*A/(4pi^2)*(2pi)*E/v^2. A = 3*sqrt(3)/2*1.42^2
    =#
    es, vs = dos(pb_model(pb_graphene.monolayer(), pb.translational_symmetry()), mesh=300, histogram_width=5)
    dosanalytic(x::Real) = 3*sqrt(3)/2*1.42^2*2*pi/(pi^2)*abs(x)/36
    for (e, v) in zip(es, vs)
        ((abs(e)<1) & (v != 0) & (dosanalytic(e) !=0)) && @test abs(100*(dosanalytic(e)-2*v)/(2*v)) < 10 
    end
end