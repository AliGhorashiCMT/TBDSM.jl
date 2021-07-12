"""
$(TYPEDSIGNATURES)
"""
function make_lattice(a1::Vector{<:Real}, a2::Vector{<:Real}; sublattices::Union{Nothing, Vector{<:Tuple{AbstractString,Vector{<:T}}}} = nothing, 
    hoppings::Union{Nothing, Vector{<:Tuple{Vector{<:Integer}, AbstractString, AbstractString, <:Real}}} = nothing ) where T
    lat = pb_lattice(a1, a2)
    if sublattices isa Vector{<:Tuple{String, Vector{T}}} where T
        num_sublattices = length(sublattices)
        for i in 1:num_sublattices
            lat.add_sublattices(sublattices[i])
        end
    end

    if hoppings isa Vector{<:Tuple{Vector{<:Integer}, AbstractString, AbstractString, <:Real}}
        for hopping in hoppings
            lat.add_hoppings(hopping)
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

