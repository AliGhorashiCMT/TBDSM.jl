"""
$(TYPEDSIGNATURES)
"""
function impol_3d(qx::Real, qy::Real, qz::Real, μ::Real, lat::PyCall.PyObject; histogram_width::Real = 10, mesh::Integer=10)
    qx = qx*10
    qy = qy*10 ##The user is expected to give wavevectors in inverse angstrom so this converts to inverse nanometers
    qz = qz*10 
    im_pols = zeros(histogram_width*30)
    b1, b2, b3 = lat.reciprocal_vectors()
    bzone_volume = abs(dot(b3, cross(b1, b2)))/1000 ## Brillouin zone volume in inverse angstrom cubed
    mod2d = pb_model(lat, pb.translational_symmetry())
    solver2d = pb_solver(mod2d)
    solver2d.set_wave_vector([0, 0])
    GammaBands = solver2d.eigenvalues
    num_bands = length(GammaBands)
    for (xiter, yiter, ziter) in CartesianIndices(rand(mesh, mesh, mesh))
        kx, ky, kz = xiter/mesh*b1[1] + yiter/mesh*b2[1], xiter/mesh*b1[2] + yiter/mesh*b2[2], ziter/mesh*b1[3]+ziter/mesh*b2[3]
        solver2d.set_wave_vector([kx, ky, kz])
        Eks = solver2d.eigenvalues
        vecks= solver2d.eigenvectors
        solver2d.set_wave_vector([kx+qx, ky+qy, kz+qz])
        Ekplusqs = solver2d.eigenvalues
        veckplusqs= solver2d.eigenvectors
        for dn_band in 1:num_bands
            for up_band in 1:num_bands
                edn = Eks[dn_band]
                vdn = vecks[:, dn_band]
                eupq = Ekplusqs[up_band]
                vupq = veckplusqs[:, up_band]
                overlap = abs(sum(conj(vupq).*vdn))^2
                ω = eupq - edn
                f2 = heaviside(μ-eupq)
                f1 = heaviside(μ-edn)
                ω>0 || continue
                im_pols[round(Int, ω*histogram_width+1)] += 2/(2π)^3*π*histogram_width*(f2-f1)*overlap*(1/mesh)^3*bzone_volume
            end
        end
    end
    return im_pols
end

"""
$(TYPEDSIGNATURES)
"""
function realeps_3d(qx::Real, qy::Real, qz::Real, ωs::Real, μ::Real, lat::PyCall.PyObject; histogram_width::Real=10, mesh::Integer=100, max_energy_integration::Real=20)
    q=abs(sqrt(qx^2+qy^2+qz^2))
    im_pol = impol_3d(qx, qy, qz, μ, lat, mesh=mesh, histogram_width=histogram_width)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(30-1/histogram_width), im_pol)
    epsilons = zeros(length(ωs))
    ErrorAbs=1e-20
    for (index, ω) in enumerate(ωs)
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
        epsilons[index] = 1-90.5*2/q^2*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
    end
    return epsilons
    
end