"""
$(TYPEDSIGNATURES)

Returns the densities of states for a given Pybinding model
"""
function dos(model::PyCall.PyObject; mesh=100, histogram_width=100)
    d = length(model.lattice.reciprocal_vectors())
    println(d)
    return dos(model, Val(d), mesh=mesh, histogram_width=histogram_width)
end

function dos(model::PyCall.PyObject, ::Val{1}; mesh=100, histogram_width=100)
    println("Dimension: ", 1)
    solver = pb_solver(model)
    reciprocal_vector = model.lattice.reciprocal_vectors()[1]
    all_energies = Float64[]
    for i in 1:mesh
        solver.set_wave_vector(reciprocal_vector*i/mesh)
        energies = solver.eigenvalues
        for energy in energies
            push!(all_energies, energy) 
        end
    end
    mine, maxe = minimum(all_energies), maximum(all_energies)
    DOS_ARRAY = zeros(Int(round(maxe-mine)*histogram_width)+1)
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh
    end
    return DOS_ARRAY
end

function dos(model::PyCall.PyObject, ::Val{2}; mesh=100, histogram_width=100)
    println("Dimension: ", 2)
    solver = pb_solver(model)
    b1, b2 = model.lattice.reciprocal_vectors()
    all_energies = Float64[]
    for i in 1:mesh
        for j in 1:mesh
            solver.set_wave_vector(b1*i/mesh+b2*j/mesh)
            energies = solver.eigenvalues
            for energy in energies
                push!(all_energies, energy) 
            end
        end
    end
    mine, maxe = minimum(all_energies), maximum(all_energies)
    DOS_ARRAY = zeros(Int(round(maxe-mine)*histogram_width)+1)
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh^2
    end
    return DOS_ARRAY
end

function dos(model::PyCall.PyObject, ::Val{3}; mesh=100, histogram_width=100)
    println("Dimension: ", 3)
    solver = pb_solver(model)
    b1, b2, b3 = model.lattice.reciprocal_vectors()
    all_energies = Float64[]
    for i in 1:mesh
        for j in 1:mesh
            for k in 1:mesh
                solver.set_wave_vector(b1*i/mesh+b2*j/mesh+b3*k/mesh)
                energies = solver.eigenvalues
                for energy in energies
                    push!(all_energies, energy) 
                end
            end
        end
    end
    mine, maxe = minimum(all_energies), maximum(all_energies)
    DOS_ARRAY = zeros(Int(round(maxe-mine)*histogram_width)+1)
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh^3
    end
    return DOS_ARRAY
end