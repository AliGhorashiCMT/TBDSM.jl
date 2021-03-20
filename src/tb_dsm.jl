module tb_dsm
#Dependencies
using PyCall, PyPlot

const pb = PyNULL()
const pb_repo = PyNULL()
const pb_graphene = PyNULL()
const model = PyNULL()
const solver = PyNULL()
export pb, pb_repo, pb_graphene, model, solver

function __init__()
    copy!(pb, pyimport("pybinding"))
    copy!(pb_repo, pyimport("pybinding.repository"))
    copy!(pb_graphene, pyimport("pybinding.repository.graphene"))
    copy!(model, pb.Model)
    copy!(solver, pb.solver.lapack)
end

include("tb_model.jl")

include("tb_2d_dielectric.jl")

include("tb_3d_dielectric.jl")

include("tb_specific_examples.jl")
export graphene_bands

end # module
