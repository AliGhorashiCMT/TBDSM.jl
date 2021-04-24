"""
$(TYPEDSIGNATURES)
"""
function plot_density(model::PyCall.PyObject; orbitals::Integer=1, mesh=10)
    #positions = [val.position[1:2] for (key, val) in model.lattice.sublattices]
    positions = [val.position[1:2] for (key, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    println(positions)
    solver = pb_solver(model)
    solver.set_wave_vector([1/2, 1/9])
    wfnssquared = abs.(solver.eigenvectors[:, orbitals]).^2
    @assert length(wfnssquared) == length(positions)
    density_array = zeros(mesh, mesh)
    for i in 1:mesh
        for j in 1:mesh
            x, y = (i-mesh/2)/mesh*4, (j-mesh/2)/mesh*4
            #for (index, pos) in enumerate(positions)
            for (roverlap, pos) in zip(wfnssquared, positions)
                x0, y0 = pos
                density_array[i, j] += exp(-(x-x0)^2*50-(y-y0)^2*50)*roverlap#wfnssquared[index]
            end
        end
    end
    println((argmax(density_array)[1]-mesh/2)/mesh*4,  "   ", (argmax(density_array)[2]-mesh/2)/mesh*4)
    display(Plots.heatmap(transpose(density_array)))
end