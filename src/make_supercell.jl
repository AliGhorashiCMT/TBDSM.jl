"""
$(TYPEDSIGNATURES)
"""
function make_supercell(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}}, cell_mult::Vector{<:Integer}, d::Real )
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    for (x, y) in Tuple.(CartesianIndices(rand(cell_mult...)))
        for (index, sublattice) in enumerate(sublattices)
            try
            lat.add_sublattices(("$(sublattice[1])$(index)$(x-1)$(y-1)", sum(lat_vectors.*[x, y])+sublattice[2] ))
            catch
                println("exception")
            end
        end
    end
    for (index1, sublattice1) in enumerate(sublattices)
        for (x1, y1) in Tuple.(CartesianIndices(rand(cell_mult...)))
            for  (index2, sublattice2) in enumerate(sublattices)
                sublatdistance = sublattice2[2] - sublattice1[2]
                for (x2, y2) in Tuple.(CartesianIndices(rand(cell_mult...)))
                    for (supercellx, supercelly) in Tuple.(CartesianIndices(rand(2, 2)))
                        try    
                            distsquared = sum((lat_vectors[1]*(x2-x1) + lat_vectors[2]*(y2-y1) + newvecs[1]*(supercellx-1) + newvecs[2]*(supercelly-1) + sublatdistance).^2)
                            distsquared ≈ d^2 && lat.add_hoppings(([supercellx-1, supercelly-1], "$(sublattice1[1])$(index1)$(x1-1)$(y1-1)", "$(sublattice2[1])$(index2)$(x2-1)$(y2-1)", 1)) 
                        catch 
                        end
                    end
                end
            end
        end
    end
    return pb_model(lat, pb.translational_symmetry()), lat
end

"""
$(TYPEDSIGNATURES)
Faster version of make_supercell- explicitly finds supercell hoppings instead of relying on nearest-neighbor distances
"""
function make_supercell2(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}}, 
    cell_mult::Vector{<:Integer}, hops::Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer}}} )
    plt.close()
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    for (x1, y1) in Tuple.(CartesianIndices(rand(cell_mult...)))
        x, y = x1-1, y1-1
        for (_, sublattice) in enumerate(sublattices)
            try
                lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2]))
            catch
                println("exception")
            end
        end
    end
    totalnumhops = 0
    for (index1, sublattice1) in enumerate(sublattices)
        for (x2, y2) in Tuple.(CartesianIndices(rand(cell_mult...)))
            x1, y1 = x2-1, y2-1
            for hop in hops
                hop[1] == index1 || continue
                a, b = return_superhoppings([x1, y1], hop[3], cell_mult)
                secondsub = sublattices[hop[2]][1]
                try
                    lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", 1))
                    totalnumhops +=1
                catch
                    println("some error")
                end
            end
        end
    end
    println("Total number of hoppings (check): ", totalnumhops)
    supermodel = pb_model(lat, pb.translational_symmetry())
    supermodel.plot(site=Dict("radius" =>0.1))
    plt.savefig("lat.png")
    return supermodel
end

function make_defectcell(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}},
    sublatonsite::Vector{<:Real}, cell_mult::Vector{<:Integer}, hops::Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer}, <:Integer}}, 
    defectenergy::Real, defecthop::Real )

    @assert length(sublattices) == length(sublatonsite)
    plt.close()
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    for (x1, y1) in Tuple.(CartesianIndices(rand(cell_mult[1], cell_mult[2])))
        x, y = x1-1, y1-1
        for (index, sublattice) in enumerate(sublattices)
            try
            (x==0 && y==0 && index==1) ? lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2], defectenergy )) : lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2], sublatonsite[index] ))
            catch e
                println(e)
                println("exception")
            end
        end
    end
    totalnumhops = 0
    for (index1, sublattice1) in enumerate(sublattices)
        for (x2, y2) in Tuple.(CartesianIndices(rand(cell_mult[1], cell_mult[2])))
            x1, y1 = x2-1, y2-1
            for hop in hops
                hop[1] == index1 || continue
                a, b = return_superhoppings([x1, y1], hop[3], cell_mult)
                secondsub = sublattices[hop[2]][1]
                try
                    ("$(sublattice1[1])$(x1)$(y1)" == sublattices[1][1]*"00" ||  "$(secondsub)$(a...)" == sublattices[1][1]*"00") ? println("Making defect hopping\n") : nothing
                    ("$(sublattice1[1])$(x1)$(y1)" == sublattices[1][1]*"00" ||  "$(secondsub)$(a...)" == sublattices[1][1]*"00") ? lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", defecthop)) : lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", hop[4]))
                    totalnumhops +=1
                catch
                    println("some error")
                end
            end
        end
    end
    println("Total number of hoppings (check): ", totalnumhops)
    supermodel = pb_model(lat, pb.translational_symmetry())
    return supermodel
end

function make_hubbarddefectcell(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}}, sublatonsite::Vector{<:Real}, cell_mult::Vector{<:Integer}, 
    hops::Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer}, <:Integer}}, defectenergy::Real, defecthop::Real, U::Real)
    @assert length(sublattices) == length(sublatonsite)
    plt.close()
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    for (x1, y1) in Tuple.(CartesianIndices(rand(cell_mult[1], cell_mult[2])))
        x, y = x1-1, y1-1
        for (index, sublattice) in enumerate(sublattices)
            try
            (x==0 && y==0 && index==1) ? lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2], [defectenergy+U, defectenergy] )) : lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2], [sublatonsite[index], sublatonsite[index]] ))
            catch e
                println(e)
                println("exception")
            end
        end
    end
    totalnumhops = 0
    for (index1, sublattice1) in enumerate(sublattices)
        for (x2, y2) in Tuple.(CartesianIndices(rand(cell_mult[1], cell_mult[2])))
            x1, y1 = x2-1, y2-1
            for hop in hops
                hop[1] == index1 || continue
                a, b = return_superhoppings([x1, y1], hop[3], cell_mult)
                secondsub = sublattices[hop[2]][1]
                try
                    ("$(sublattice1[1])$(x1)$(y1)" == sublattices[1][1]*"00" ||  "$(secondsub)$(a...)" == sublattices[1][1]*"00") && println("Making defect hopping\n") 
                    ("$(sublattice1[1])$(x1)$(y1)" == sublattices[1][1]*"00" ||  "$(secondsub)$(a...)" == sublattices[1][1]*"00") ? lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", [[defecthop, 0], [0, defecthop]])) : lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", [[hop[4], 0], [0, hop[4]]]))
                    totalnumhops +=1
                catch
                    println("some error")
                end
            end
        end
    end
    println("Total number of hoppings (check): ", totalnumhops)
    supermodel = pb_model(lat, pb.translational_symmetry())
    return supermodel
end


"""
$(TYPEDSIGNATURES)
Returns many defect cells using make_hubbarddefectcell. The idea is to be able to make many different 
defect cells to fit to Ab Initio data.
"""
function make_hubbarddefectcells(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}}, sublatonsites::Vector{<:Vector{<:Real}}, cell_mult::Vector{<:Integer}, 
    hopses::Vector{<:Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer}, <:Integer}}}, defectenergies::Vector{<:Real}, defecthops::Vector{<:Real}, Us::Vector{<:Real})
    defectcells = PyCall.PyObject[] 
    for sublatonsite in sublatonsites
        for hops in hopses
            for defectenergy in defectenergies
                for defecthop in defecthops
                    for U in Us
                        push!(defectcells, make_hubbarddefectcell(lat_vectors, sublattices, sublatonsite, cell_mult, hops, defectenergy, defecthop, U))
                    end
                end
            end
        end
    end
    return defectcells
end

"""
$(TYPEDSIGNATURES)
"""
function make_hubbarddefectcells(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:AbstractString, Vector{<:Real}}}, sublatonsite::Vector{<:Real}, cell_mult::Vector{<:Integer}, 
    hops::Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer}, <:Integer}}, defectenergies::Vector{<:Real}, defecthops::Vector{<:Real}, Us::Vector{<:Real})

    make_hubbarddefectcells(lat_vectors, sublattices, [sublatonsite], cell_mult, 
        [hops], defectenergies, defecthops, Us)
end
"""
$(TYPEDSIGNATURES)
"""
function make_onedsupercell(supercellmult::Integer)
    N = supercellmult
    println("Supercell bands for multiplicity: ", N)
    oned_lat = pb_lattice([N])
    oned_lat.add_sublattices([("A$(n)", n) for n in 0:N-1]...) #Add all atoms within unit cell
    oned_lat.add_hoppings([([0], "A$(n-1)", "A$(n)", 1) for n in 1:N-1]...) #Add all intra unit cell hoppings
    oned_lat.add_hoppings(([-1], "A0", "A$(N-1)", 1)) #Add the one inter unit cell hopping
    oned_mod = pb_model(oned_lat, pb.translational_symmetry())
    Energies = Array{Float64, 2}(undef, (length(-π/N:π/(100*N):π/N), N))
    Eigenvectors = Array{ComplexF64, 3}(undef, (length(-π/N:π/(100*N):π/N), N, N)) #N component eigenvectors and N bands 
    for (kindex, x) in enumerate(-π/N:π/(100*N):π/N) #Only look at first Brillouin Zone
        oned_solver = pb_solver(oned_mod)
        oned_solver.set_wave_vector(x)
        Energies[kindex, :] = oned_solver.eigenvalues
        Eigenvectors[kindex, :, :] = oned_solver.eigenvectors
    end
    return Energies, Eigenvectors
end

"""
$(TYPEDSIGNATURES)
"""
function make_supercellbands(supercellmult::Integer)
    N = supercellmult
    println("Supercell bands for multiplicity: ", N)
    oned_lat = pb_lattice([1])
    oned_lat.add_sublattices(("A", [0, 0]))
    oned_lat.add_hoppings(([1], "A", "A", 1))
    oned_mod = pb_model(oned_lat, pb.translational_symmetry())
    Energies = Array{Float64, 2}(undef, (length(-π/N:π/(100*N):π/N), N))
    for (kindex, x) in enumerate(-π/N:π/(100*N):π/N)
        for n in 1:N
            oned_solver = pb_solver(oned_mod)
            oned_solver.set_wave_vector(x+(n-1)*2π/N)
            Energies[kindex, n] = oned_solver.eigenvalues[1]
        end
    end
    return Energies
end

function make_graphenesupercellbands(supercellmult::Integer)
    N = supercellmult
    println("Supercell bands for multiplicity: ", N)

    lat = pb_lattice([1, 0], [1/2, 1/2*sqrt(3)])
    lat.add_sublattices(
        ("A", [0, -1/(2*sqrt(3))]),
        ("B", [0,  1/(2*sqrt(3))])
    )
    lat.add_hoppings(
        ([0,  0], "A", "B", 1),
        ([1, -1], "A", "B", 1),
        ([0, -1], "A", "B", 1))
    
    K2 = [2*pi / (3), 2*pi / sqrt(3)]
    M = [0, 2π/sqrt(3)]
    K1 = [-4π/3, 0]
    kpath = pb.make_path(K1, [0, 0], M, K2, [0, 0], -K1,  step=0.02)
    graphene_mod = pb_model(lat, pb.translational_symmetry())
    Energies = Array{Float64, 2}(undef, (size(kpath)[1], 2N^2))
    b1, b2 = lat.reciprocal_vectors()
    b1 *= 1/N
    b2 *= 1/N
    println(b1, b2)
    for (kindex, x) in enumerate(eachrow(kpath))
        k = collect(x)
        allsups = Array{Float64, 3}(undef, (N, N, 2))
        for nx in 0:N-1
            for ny in 0:N-1
                solver = pb_solver(graphene_mod)
                solver.set_wave_vector(k.+nx*b1[1:2].+ny*b2[1:2])
                allsups[nx+1, ny+1, 1] = solver.eigenvalues[1]
                allsups[nx+1, ny+1, 2] = solver.eigenvalues[2]
            end
        end
        Energies[kindex, :] = reshape(allsups, 2N^2)
    end
    return Energies
end

function return_superhoppings(superid::Vector{<:Integer}, hopping::Vector{<:Integer}, supercell::Vector{<:Integer})
    supercellhopping = div.((hopping+superid), supercell, RoundDown)
    superid2 = rem.((hopping+superid), supercell)
    superid2[1]<0 ? superid2[1] = supercell[1]+ superid2[1] : nothing
    superid2[2]<0 ? superid2[2] = supercell[2]+ superid2[2] : nothing
    return superid2, supercellhopping
end