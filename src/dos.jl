"""
$(TYPEDSIGNATURES)

Returns the densities of states for a given Pybinding model
"""
function dos(model::PyCall.PyObject; mesh=100, histogram_width=100)
    d = length(model.lattice.reciprocal_vectors())
    println(d)
    return dos(model, Val(d), mesh=mesh, histogram_width=histogram_width)
end

"""
$(TYPEDSIGNATURES)

"""
function bands_overlayeddos(model::PyCall.PyObject, ks...; mesh=100, histogram_width=100, normalized::Bool=true)
    d = length(model.lattice.reciprocal_vectors())
    println(d)
    return bands_overlayeddos(model, Val(d), ks..., mesh=mesh, histogram_width=histogram_width, normalized=normalized)
end

function bands_overlayeddos(model::PyCall.PyObject, ::Val{1},  ks...; normalized::Bool=true, mesh=100, histogram_width=100)
    solver = pb_solver(model)
    diffks = diff(float.(collect(ks)))
    bands = Vector{Vector{Float64}}()
    for (k, diffk) in zip(ks, diffks)
        for n in 1:100
            kinterpolate = k + diffk*n/100
            solver.set_wave_vector(kinterpolate)
            push!(bands, solver.eigenvalues)
        end
    end
    bands_2d = zeros(length(bands), length(bands[1]) )
    for row in 1:length(bands)
        bands_2d[row, :]=bands[row]
    end
    println(bands)
    println(bands_2d)
    a = Plots.plot(bands_2d, legend=false)
    dose, dosv = dos(model, mesh=mesh)
    b = Plots.plot(dosv, dose, legend=false)
    display(Plots.plot(a, b, size=(1000, 500), linewidth=5))
end

function bands_overlayeddos(model::PyCall.PyObject, ::Val{2},  ks...; normalized::Bool=true, mesh=100, histogram_width=100, kwargs...)
    solver = pb_solver(model)
    diffks = diff(float.(collect(ks)))
    bands = Vector{Vector{Float64}}()
    b1, b2 = model.lattice.reciprocal_vectors()
    for (k, diffk) in zip(ks, diffks)
        for n in 1:100
            kinterpolate = k + diffk*n/100
            normalized ? solver.set_wave_vector(sum(kinterpolate.*[b1, b2])) : solver.set_wave_vector(kinterpolate)
            push!(bands, solver.eigenvalues)
        end
    end
    bands_2d = zeros(length(bands), length(bands[1]) )
    println(length(bands[1]))
    for row in 1:length(bands)
        bands_2d[row, :]=bands[row]
    end
    a = Plots.plot(bands_2d, legend=false)
    dose, dosv = dos(model, mesh=mesh, histogram_width=histogram_width)
    b = Plots.plot(dosv, dose, legend=false)
    display(Plots.plot(a, b, size=(1000, 500), linewidth=5, xticks=[]))
end

function dos(model::PyCall.PyObject, ::Val{1}; nelectrons::Union{Nothing, Real}=nothing, mesh=100, histogram_width=100)
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
    totlength= Int(round(maxe-mine)*histogram_width)+1
    DOS_ARRAY = zeros(totlength)
    E_ARRAY = mine:(maxe-mine)/totlength:maxe-(maxe-mine)/totlength#collect(mine:1/histogram_width:(length(DOS_ARRAY)-1)/histogram_width+mine)
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh
    end
    if nelectrons isa Real  
        println(sum(diff(E_ARRAY).*DOS_ARRAY[1:end-1]))
        @assert sum(diff(E_ARRAY).*DOS_ARRAY[1:end-1]) ≈ nelectrons
    end
    @assert length(E_ARRAY) == length(DOS_ARRAY)
    return E_ARRAY, DOS_ARRAY
end

function dos(model::PyCall.PyObject, ::Val{2}; mesh=100, histogram_width=100)
    println("Dimension: ", 2)
    solver = pb_solver(model)
    b1, b2 = model.lattice.reciprocal_vectors()
    all_energies = Float64[]
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        solver.set_wave_vector(b1*i/mesh+b2*j/mesh)
        energies = solver.eigenvalues
        for energy in energies
            push!(all_energies, energy) 
        end
    end
    mine, maxe = minimum(all_energies), maximum(all_energies)
    println(mine, " ", maxe)
    totlength = Int(round((maxe-mine)*histogram_width))+1
    E_ARRAY = range(mine, maxe, length=totlength)
    DOS_ARRAY = zeros(totlength)
    @assert length(E_ARRAY) == length(DOS_ARRAY) 
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh^2
    end
    #TODO @assert sum(diff(E_ARRAY).*DOS_ARRAY[1:end-1])
    return E_ARRAY, DOS_ARRAY
end

function dos(model::PyCall.PyObject, ::Val{3}; mesh=100, histogram_width=100)
    println("Dimension: ", 3)
    solver = pb_solver(model)
    b1, b2, b3 = model.lattice.reciprocal_vectors()
    all_energies = Float64[]
    for (i, j, k) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        solver.set_wave_vector(b1*i/mesh+b2*j/mesh+b3*k/mesh)
        energies = solver.eigenvalues
        for energy in energies
            push!(all_energies, energy) 
        end
    end
    mine, maxe = minimum(all_energies), maximum(all_energies)
    DOS_ARRAY = zeros(Int(round(maxe-mine)*histogram_width)+1)
    for energy in all_energies
        DOS_ARRAY[Int(round((energy-mine)*histogram_width))+1] += histogram_width/mesh^3
    end
    return DOS_ARRAY
end