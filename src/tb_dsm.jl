module tb_dsm
#Dependencies
using PyCall, PyPlot

const pb = PyNULL()
const pb_repo = PyNULL()
const pb_graphene = PyNULL()
const pb_model = PyNULL()
const pb_solver = PyNULL()
const pb_lattice = PyNULL()
export pb, pb_repo, pb_graphene, pb_model, pb_solver, pb_lattice

function __init__()
    copy!(pb, pyimport("pybinding"))
    copy!(pb_repo, pyimport("pybinding.repository"))
    copy!(pb_graphene, pyimport("pybinding.repository.graphene"))
    copy!(pb_model, pb.Model)
    copy!(pb_solver, pb.solver.lapack)
    copy!(pb_lattice, pb.Lattice)
end

include("tb_model.jl")

include("tb_2d_dielectric.jl")

include("tb_3d_dielectric.jl")

include("tb_specific_examples.jl")
export graphene_bands

end # module
