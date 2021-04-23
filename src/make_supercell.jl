"""
$(TYPEDSIGNATURES)
"""
function make_supercell(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:String, Vector{<:Real}}}, cell_mult::Vector{<:Integer}, d::Real )
    plt.close()
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    #lat.add_sublattices(("A", [0, 0]))
    for x in 0:cell_mult[1]-1
        for y in 0:cell_mult[2]-1
            for (index, sublattice) in enumerate(sublattices)
                try
                lat.add_sublattices(("$(sublattice[1])$(index)$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2] ))
                catch
                    println("exception")
                end
            end
        end
    end
    for (index1, sublattice1) in enumerate(sublattices)
    for x1 in 0:cell_mult[1]-1
        for y1 in 0:cell_mult[2]-1
            for  (index2, sublattice2) in enumerate(sublattices)
            sublatdistance = sublattice2[2] - sublattice1[2]
            for x2 in 0:cell_mult[1]-1
                for y2 in 0:cell_mult[2]-1  
                    for supercellx in -1:1
                        for supercelly in -1:1  
                            try    
                                distsquared = sum((lat_vectors[1]*(x2-x1) + lat_vectors[2]*(y2-y1) + newvecs[1]*supercellx + newvecs[2]*supercelly + sublatdistance).^2)
                                distsquared ≈ d^2 ? lat.add_hoppings(([supercellx, supercelly], "$(sublattice1[1])$(index1)$(x1)$(y1)", "$(sublattice2[1])$(index2)$(x2)$(y2)", 1)) : nothing
                            catch e
                            end
                        end
                    end
                end
            end
        end
    end
    end
    end
    #lat.plot()
    return pb_model(lat, pb.translational_symmetry())
    #pb_solver(graphene_mod).calc_bands(K1, Gamma, M, K2).plot()
   # plt.show()
end

"""
$(TYPEDSIGNATURES)
Faster version of make_supercell- explicitly finds supercell hoppings instead of relying on nearest-neighbor distances
"""
function make_supercell2(lat_vectors::Vector{<:Vector{<:Real}}, sublattices::Vector{<:Tuple{<:String, Vector{<:Real}}}, cell_mult::Vector{<:Integer}, hops::Vector{<:Tuple{<:Integer, <:Integer, <:Vector{<:Integer} }} )
    plt.close()
    newvecs = lat_vectors.*cell_mult
    println(newvecs[1])
    lat = pb_lattice(newvecs...) # pybinding lattice with superlattice vectors
    for x in 0:cell_mult[1]-1
        for y in 0:cell_mult[2]-1
            for (index, sublattice) in enumerate(sublattices)
                try
                lat.add_sublattices(("$(sublattice[1])$(x)$(y)", sum(lat_vectors.*[x, y])+sublattice[2] ))
                catch
                    println("exception")
                end
            end
        end
    end
    totalnumhops = 0
    for (index1, sublattice1) in enumerate(sublattices)
        for x1 in 0:cell_mult[1]-1
            for y1 in 0:cell_mult[2]-1
                for hop in hops
                    if hop[1] == index1 
                        a, b = return_superhoppings([x1, y1], hop[3], cell_mult)
                        secondsub = sublattices[hop[2]][1]
                        try
                            lat.add_hoppings((b, "$(sublattice1[1])$(x1)$(y1)", "$(secondsub)$(a...)", 1))
                            totalnumhops +=1
                        catch
                            println("some error")
                            @warn "some error"
                        end
                    else
                        continue
                    end
                end
            end
        end
    end
    println("Total number of hoppings (check): ", totalnumhops)
    supermodel = pb_model(lat, pb.translational_symmetry())
    #supermodel.plot(site=Dict("radius" =>1))
    #plt.savefig("lat.png")
    #plt.close()
    return supermodel
end

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