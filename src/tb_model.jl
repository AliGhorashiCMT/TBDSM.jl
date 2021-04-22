"""
$(TYPEDSIGNATURES)
"""
function make_lattice(a1::Array{<:Real, 1}, a2::Array{<:Real, 1}; sublattices::Union{Nothing, Array{Tuple{String,Array{T,1}},1}} = nothing, hoppings::Union{Nothing, Array{<:Tuple{Array{Int64,1},String,String, <:Real},1}} = nothing ) where T
    lat = pb_lattice(a1, a2)
    if sublattices isa Array{Tuple{String,Array{T,1}},1} where T
        num_sublattices = length(sublattices)
        for i in 1:num_sublattices
            lat.add_sublattices(sublattices[i])
        end
    end

    if hoppings isa Array{<:Tuple{Array{Int64,1},String,String, <:Real},1}
        num_hoppings = length(hoppings)
        for i in 1:num_hoppings
            lat.add_hoppings(hoppings[i])
        end
    end
    return lat
end

"""
$(TYPEDSIGNATURES)
"""
function plot_bands(lat::PyCall.PyObject, kvecs...)
    model = pb_model(lat, pb.translational_symmetry())
    solver = pb_solver(model)
    return solver.calc_bands(kvecs...) ##Will return the band structure along the specified kvector path
end

