"""
$(TYPEDSIGNATURES)
"""
function plot_wfn(model::PyCall.PyObject; orbital::Integer=1, mesh=10)
    #positions = [val.position[1:2] for (key, val) in model.lattice.sublattices]
    positions = [val.position[1:2] for (key, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    println(positions)
    solver = pb_solver(model)
    solver.set_wave_vector([1/2, 1/9])
    wfnssquared = abs.(solver.eigenvectors[:, orbital]).^2
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

"""
$(TYPEDSIGNATURES)
Returns a density heatmap of the supplied Pybinding model (for electrons up to the supplied Fermi energy (μ)). Orbitals are modeled as normalized 2d Gaussians.
"""
function plot_density(model::PyCall.PyObject; μ::Real=100, mesh::Integer=100, bzonemesh::Integer=20, σ=1/10, extent::Integer=5)
    positions = [val.position[1:2] for (key, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    println(positions)
    solver = pb_solver(model)
    density_vec = zeros(length(positions))
    for i in 1:bzonemesh
        for j in 1:bzonemesh
            solver.set_wave_vector([i/bzonemesh, j/bzonemesh])
            for (energy, orbital) in zip(solver.eigenvalues, eachcol(solver.eigenvectors))
                density_vec += (heaviside(μ-energy)/(bzonemesh^2))*(abs.(orbital)).^2
            end
        end
    end
    density_array = zeros(mesh, mesh)
    for i in 1:mesh
        for j in 1:mesh
            x, y = (i-mesh/2)/mesh*extent, (j-mesh/2)/mesh*extent
            for (density, pos) in zip(density_vec, positions)
                x0, y0 = pos
                density_array[i, j] += 1/(2*π*σ^2)*exp(-((x-x0)^2+(y-y0)^2)/(2*σ^2))*density
            end
        end
    end
    println((argmax(density_array)[1]-mesh/2)/mesh*extent,  "   ", (argmax(density_array)[2]-mesh/2)/mesh*extent)
    println("Number of electrons in density plot: ", sum(density_array*(extent/mesh)^2))
    display(Plots.heatmap(transpose(density_array)))
end