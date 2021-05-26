"""
$(TYPEDSIGNATURES)
"""
function plot_wfn(model::PyCall.PyObject; orbital::Integer=1, mesh::Integer=10)
    positions = [val.position[1:2] for (_, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    println(positions)
    solver = pb_solver(model)
    solver.set_wave_vector([1/2, 1/9])
    wfnssquared = abs.(solver.eigenvectors[:, orbital]).^2
    @assert length(wfnssquared) == length(positions)
    density_array = zeros(mesh, mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(bzonemesh, bzonemesh)))
        x, y = (i-mesh/2)/mesh*4, (j-mesh/2)/mesh*4
            #for (index, pos) in enumerate(positions)
        for (roverlap, pos) in zip(wfnssquared, positions)
            x0, y0 = pos
            density_array[i, j] += exp(-(x-x0)^2*50-(y-y0)^2*50)*roverlap#wfnssquared[index]
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
    positions = [val.position[1:2] for (_, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    println(positions)
    solver = pb_solver(model)
    density_vec = zeros(length(positions))
    for (i, j) in Tuple.(CartesianIndices(rand(bzonemesh, bzonemesh)))
        solver.set_wave_vector([i/bzonemesh, j/bzonemesh])
        for (energy, orbital) in zip(solver.eigenvalues, eachcol(solver.eigenvectors))
            density_vec += (heaviside(μ-energy)/(bzonemesh^2))*(abs.(orbital)).^2
        end
    end
    density_array = zeros(mesh, mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        x, y = (i-mesh/2)/mesh*extent, (j-mesh/2)/mesh*extent
        for (density, pos) in zip(density_vec, positions)
            x0, y0 = pos
            density_array[i, j] += 1/(2*π*σ^2)*exp(-((x-x0)^2+(y-y0)^2)/(2*σ^2))*density
        end
    end
    println((argmax(density_array)[1]-mesh/2)/mesh*extent,  "   ", (argmax(density_array)[2]-mesh/2)/mesh*extent)
    println("Number of electrons in density plot: ", sum(density_array*(extent/mesh)^2))
    display(Plots.heatmap(transpose(density_array)))
end

"""
$(TYPEDSIGNATURES)

Plots the spin polarized density of a given Pybinding model
"""
function plot_spinpolarizeddensity(model::PyCall.PyObject; μ::Real=100, mesh::Integer=100, bzonemesh::Integer=20, σ=1/10, extent::Integer=5)
    positions1 = [val.position[1:2] for (_, val) in sort(collect(model.lattice.sublattices), by=x->x[1])] #Need to sort dictionary by order in which sublattices were added
    positions = similar(positions1, 0)
    for pos in positions1
        push!(positions, pos, pos)
    end
    println(positions)
    solver = pb_solver(model)
    density_vecup = zeros(length(positions))
    density_vecdn = zeros(length(positions))
    for (i, j) in Tuple.(CartesianIndices(rand(bzonemesh, bzonemesh)))
        solver.set_wave_vector([i/bzonemesh, j/bzonemesh])
        numbands = length(solver.eigenvalues)
        evenbands = numbands%2==0 ? collect(2:2:numbands) :  collect(2:2:numbands-1)
        oddbands = numbands%2==0 ? collect(1:2:numbands-1) :  collect(1:2:numbands)

        for (energy, orbital) in zip(solver.eigenvalues, eachcol(solver.eigenvectors))
            sum(round.(collect(orbital)[evenbands], digits=2) .== 0) == 0 && (density_vecup += (heaviside(μ-energy)/(bzonemesh^2))*(abs.(orbital)).^2)
            sum(round.(collect(orbital)[oddbands], digits=2) .== 0) == 0 && (density_vecdn += (heaviside(μ-energy)/(bzonemesh^2))*(abs.(orbital)).^2)            
        end
    end
    println("up: ", sum(density_vecup))
    println("dn: ", sum(density_vecdn))

    println(length(density_vecdn))
    
    density_arrayup = zeros(mesh, mesh)
    density_arraydn = zeros(mesh, mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        x, y = (i-mesh/2)/mesh*extent, (j-mesh/2)/mesh*extent
        for (densityup, densitydn, pos) in zip(density_vecup, density_vecdn, positions)
            x0, y0 = pos
            density_arrayup[i, j] += 1/(2*π*σ^2)*exp(-((x-x0)^2+(y-y0)^2)/(2*σ^2))*(densityup)
            density_arraydn[i, j] += 1/(2*π*σ^2)*exp(-((x-x0)^2+(y-y0)^2)/(2*σ^2))*(densitydn)
        end
    end
    display(Plots.heatmap(transpose(density_arraydn-density_arrayup)))
end