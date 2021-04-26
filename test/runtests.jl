using Pkg;
ENV["PYTHON"]=""
Pkg.build("PyCall")
using Test, tb_dsm

@testset "First test" begin
    @test 1==1
end