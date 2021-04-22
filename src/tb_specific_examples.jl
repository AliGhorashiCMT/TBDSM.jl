"""
$(TYPEDSIGNATURES)
"""
function graphene_bands()
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    graphene_mod = pb_model(pb_graphene.monolayer(), pb.translational_symmetry())
    pb_solver(graphene_mod).calc_bands(K1, Gamma, M, K2).plot()
end

"""
$(TYPEDSIGNATURES)
"""
function heaviside(x::Real)
    x>0 ? 1 : 0
end

"""
$(TYPEDSIGNATURES)
"""
function graphene_impol(qx::Real, qy::Real, μ::Real;mesh::Int=10, histogram_width::Real=100)
    qx = qx*10
    qy = qy*10 ##The user is expected to give wavevectors in inverse angstrom so this converts to inverse nanometers
    im_pols = zeros(histogram_width*30)
    b1, b2 = pb_graphene.monolayer().reciprocal_vectors()
    bzone_area = abs(cross(b1, b2)[3])/100 ## Brillouin zone area in inverse angstrom squared
    graphene_mod = pb_model(pb_graphene.monolayer(), pb.translational_symmetry())
    graphene_solver = pb_solver(graphene_mod)
    for xiter in 1:mesh
        for yiter in 1:mesh
            kx, ky = xiter/mesh*b1[1] + yiter/mesh*b2[1], xiter/mesh*b1[2] + yiter/mesh*b2[2]
            graphene_solver.set_wave_vector([kx, ky])
            Ednk, Eupk = graphene_solver.eigenvalues[1], graphene_solver.eigenvalues[2]
            vecdnk, vecupk = graphene_solver.eigenvectors[:, 1], graphene_solver.eigenvectors[:, 2]
            graphene_solver.set_wave_vector([kx, ky]+[qx, qy])
            Ednkpluq, Eupkplusq = graphene_solver.eigenvalues[1], graphene_solver.eigenvalues[2]
            vecdnkplusq, vecupkplusq = graphene_solver.eigenvectors[:, 1], graphene_solver.eigenvectors[:, 2]
            overlap1 = abs(sum(conj(vecupkplusq).*vecdnk))^2
            overlap2 = abs(sum(conj(vecupkplusq).*vecupk))^2
            ω = Eupkplusq - Ednk
            f2 = heaviside(μ-Eupkplusq)
            f1 = heaviside(μ-Ednk)
            im_pols[round(Int, ω*histogram_width+1)] = im_pols[round(Int, ω*histogram_width+1)] + 2/(2π)^2*π*histogram_width*(f2-f1)*overlap1*(1/mesh)^2*bzone_area
            ω = Eupkplusq - Eupk
            f1 = heaviside(μ-Eupk)
            if ω > 0
                im_pols[round(Int, ω*histogram_width+1)] = im_pols[round(Int, ω*histogram_width+1)] + 2/(2π)^2*π*histogram_width*(f2-f1)*overlap2*(1/mesh)^2*bzone_area
            end            
        end
    end
    return im_pols
end

function levitov_tbgraphene_impol(qx::Real, qy::Real; μ::Real=1.81e-3, mesh::Int=10, histogram_width::Real=10000)
    qx = qx*10
    qy = qy*10 ##The user is expected to give wavevectors in inverse angstrom so this converts to inverse nanometers
    im_pols = zeros(histogram_width*1)
    b1, b2 = levitov_tbg_model().reciprocal_vectors()
    bzone_area = abs(cross(b1, b2)[3])/100 ## Brillouin zone area in inverse angstrom squared
    tbgraphene_mod = pb_model(levitov_tbg_model(), pb.translational_symmetry())
    tbgraphene_solver = pb_solver(tbgraphene_mod)
    for xiter in 1:mesh
        for yiter in 1:mesh
            kx, ky = xiter/mesh*b1[1] + yiter/mesh*b2[1], xiter/mesh*b1[2] + yiter/mesh*b2[2]
            tbgraphene_solver.set_wave_vector([kx, ky])
            Ednk, Eupk = tbgraphene_solver.eigenvalues[1], tbgraphene_solver.eigenvalues[2]
            vecdnk, vecupk = tbgraphene_solver.eigenvectors[:, 1], tbgraphene_solver.eigenvectors[:, 2]

            tbgraphene_solver.set_wave_vector([kx, ky]+[qx, qy])
            Ednkpluq, Eupkplusq = tbgraphene_solver.eigenvalues[1], tbgraphene_solver.eigenvalues[2]

            vecdnkplusq, vecupkplusq = tbgraphene_solver.eigenvectors[:, 1], tbgraphene_solver.eigenvectors[:, 2]

            overlap1 = abs(sum(conj(vecupkplusq).*vecdnk))^2
            overlap2 = abs(sum(conj(vecupkplusq).*vecupk))^2

            ω = Eupkplusq - Ednk
            f2 = heaviside(μ-Eupkplusq)
            f1 = heaviside(μ-Ednk)
            im_pols[round(Int, ω*histogram_width+1)] = im_pols[round(Int, ω*histogram_width+1)] + 4/(2π)^2*π*histogram_width*(f2-f1)*overlap1*(1/mesh)^2*bzone_area

            ω = Eupkplusq - Eupk
            f1 = heaviside(μ-Eupk)
            if ω > 0
                im_pols[round(Int, ω*histogram_width+1)] = im_pols[round(Int, ω*histogram_width+1)] + 4/(2π)^2*π*histogram_width*(f2-f1)*overlap2*(1/mesh)^2*bzone_area
            end            
        end
    end
    return im_pols
end

function graphene_realeps(qx::Real, qy::Real, ω::Real, μ::Real; mesh::Real=300, histogram_width::Real=10, max_energy_integration::Real=20, kwargs...) 
    q=abs(sqrt(qx^2+qy^2))
    im_pol = graphene_impol(qx, qy, μ, mesh=mesh, histogram_width=histogram_width)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(30-1/histogram_width), im_pol)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    return 1-90.5/q*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end

function tbgraphene_realeps(qx::Real, qy::Real, ω::Real; μ::Real=1.81e-3, mesh::Real=300, histogram_width::Real=10000, max_energy_integration::Real=.2, kwargs...) 
    q=abs(sqrt(qx^2+qy^2))
    im_pol = levitov_tbgraphene_impol(qx, qy, mesh=mesh, histogram_width=histogram_width)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(1-1/histogram_width), im_pol)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    return 1-90.5/(12.12*q)*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end

function read_grapheneplasmon()
    read_plasmon = readdlm(joinpath(@__DIR__, "../data/grapheneplasmon.txt"), '\t', Float64, '\n')
    return transpose(read_plasmon)
end

function read_levitov_tbgrapheneplasmon()
    read_plasmon = readdlm(joinpath(@__DIR__, "../data/levitov_tbgraphene.txt"), '\t', Float64, '\n')
    return transpose(read_plasmon)
end

function read_grapheneplasmonline()
    read_plasmon = read_grapheneplasmon()
    plasmon = []
    for i in 1:50
        push!(plasmon, argmin(read_plasmon[:, i])*2/96)
    end
    return plasmon
end

function graphene_realeps(qx::Real, qy::Real, ωs::Array{<:Real, 1}, μ::Real; mesh::Real=300, histogram_width::Real=10, max_energy_integration::Real=20, kwargs...) 
    q=abs(sqrt(qx^2+qy^2))
    im_pol = graphene_impol(qx, qy, μ, mesh=mesh, histogram_width=histogram_width)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(30-1/histogram_width), im_pol)
    epsilons = zeros(length(ωs))
    ErrorAbs=1e-20
    for (index, ω) in enumerate(ωs)
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
        epsilons[index] = 1-90.5/q*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
    end
    return epsilons
end

function tbgraphene_realeps(qx::Real, qy::Real, ωs::Array{<:Real, 1}; μ::Real=1.81e-3, mesh::Real=300, histogram_width::Real=10000, max_energy_integration::Real=.2, kwargs...) 
    q=abs(sqrt(qx^2+qy^2))
    im_pol = levitov_tbgraphene_impol(qx, qy, mesh=mesh, histogram_width=histogram_width)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(1-1/histogram_width), im_pol)
    epsilons = zeros(length(ωs))
    ErrorAbs=1e-20
    for (index, ω) in enumerate(ωs)
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
        epsilons[index] = 1-90.5/(12.12*q)*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
    end
    return epsilons
end

function bilayer_graphene_bands()
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    graphene_mod = pb_model(pb_graphene.bilayer(), pb.translational_symmetry())
    pb_solver(graphene_mod).calc_bands(K1, Gamma, M, K2).plot()
end

function tmd_mo_s2()
    mos2_mod = pb_model(pb_mos2, pb.translational_symmetry())
    kpoints1, kpoints2, kpoints3 = [0, 0], mos2_mod.lattice.brillouin_zone()[1], mos2_mod.lattice.brillouin_zone()[2]
    pb_solver(mos2_mod).calc_bands(kpoints1, kpoints2, kpoints3).plot()
end

function levitov_tbg_model()
    a = 13.4   # [nm] unit cell length
    a_cc = a/sqrt(3)  # [nm] carbon-carbon distance
    t = 3.75/3*1e-3      # [eV] nearest neighbour hopping
    lat = pb_lattice([a, 0], [a/2, a/2 * sqrt(3)])
    lat.add_sublattices(
        ("A", [0, -a_cc/2]),
        ("B", [0,  a_cc/2])
    )
    lat.add_hoppings(
        ([0,  0], "A", "B", t),
        ([1, -1], "A", "B", t),
        ([0, -1], "A", "B", t))
    return lat
end
