using Pkg;
ENV["PYTHON"]=""
Pkg.build("PyCall")
using Test, TBDSM

@testset "First test" begin
    @test 1==1
end