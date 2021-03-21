function impol_2d(qx::Real, qy::Real, μ::Real, lat::PyCall.PyObject; histogram_width::Real = 10, mesh::Int=10)
    
    qx = qx*10
    qy = qy*10 ##The user is expected to give wavevectors in inverse angstrom so this converts to inverse nanometers
    
    im_pols = zeros(histogram_width*30)
    b1, b2 = lat.reciprocal_vectors()

    bzone_area = abs(cross(b1, b2)[3])/100 ## Brillouin zone area in inverse angstrom squared
    
    mod2d = pb_model(lat, pb.translational_symmetry())
    solver2d = pb_solver(mod2d)
    
    solver2d.set_wave_vector([0, 0])
    GammaBands = solver2d.eigenvalues
    num_bands = length(GammaBands)
    #println(num_bands)
    for xiter in 1:mesh
        for yiter in 1:mesh
            kx, ky = xiter/mesh*b1[1] + yiter/mesh*b2[1], xiter/mesh*b1[2] + yiter/mesh*b2[2]
            solver2d.set_wave_vector([kx, ky])
            Eks = solver2d.eigenvalues
            vecks= solver2d.eigenvectors

            solver2d.set_wave_vector([kx+qx, ky+qy])
            Ekplusqs = solver2d.eigenvalues
            veckplusqs= solver2d.eigenvectors

            for dn_band in 1:num_bands
                for up_band in 1:num_bands
                    edn = Eks[dn_band]
                    eup = Eks[up_band]
                    vdn = vecks[:, dn_band]
                    vup = vecks[:, up_band]

                    ednq = Ekplusqs[dn_band]
                    eupq = Ekplusqs[up_band]
                    vdnq = veckplusqs[:, dn_band]
                    vupq = veckplusqs[:, up_band]

                    overlap = abs(sum(conj(vupq).*vdn))^2

                    ω = eupq - edn
                    f2 = heaviside(μ-eupq)
                    f1 = heaviside(μ-edn)
                    if ω>0
                        #println(ω)
                        im_pols[round(Int, ω*histogram_width+1)] = im_pols[round(Int, ω*histogram_width+1)] + 2/(2π)^2*π*histogram_width*(f2-f1)*overlap*(1/mesh)^2*bzone_area
                    end
                end
            end
        end
    end
    return im_pols
end

function realeps_2d(qx::Real, qy::Real, ωs::Real, μ::Real, lat::PyCall.PyObject; histogram_width::Real=10, mesh::Int=100)
    q=abs(sqrt(qx^2+qy^2))
    im_pol = impol_2d(qx, qy, μ, lat, mesh=mesh, histogram_width=histogram_width)

    interpolated_ims=interpol.interp1d(0:1/histogram_width:(30-1/histogram_width), im_pol)
    epsilons = zeros(length(ωs))
    ErrorAbs=1e-20
    for (index, ω) in enumerate(ωs)
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    
        epsilons[index] = 1-90.5/q*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
    end
    return epsilons
    
end