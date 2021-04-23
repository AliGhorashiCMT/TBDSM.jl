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