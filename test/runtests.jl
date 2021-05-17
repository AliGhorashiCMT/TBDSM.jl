using Pkg;
ENV["PYTHON"]=""
Pkg.build("PyCall")
using Test, TBDSM

@testset "Specific Examples" begin
    @test TBDSM.heaviside(1) == 1
    @test TBDSM.heaviside(-1) == 0
end