"""
$(TYPEDSIGNATURES)

Return the momentum matrix elements at the point k. The vector q determines the direction of the momentum as well as the 
accuracy with which the derivative of the Bloch vector is calculated. Giving only one band index returns the intraband momentum 
matrix element, which is a function solely of the dispersion. Giving two band indices gives the 
"""
function momentum(k::Vector{<:Real}, q::Vector{<:Real},  band1::Integer, band2::Integer, lat::PyCall.PyObject)
    prefactor = (0.5*1e6/(3e18)^2)/(6.6e-16)
    kx, ky = k
    qx, qy = q

    mod2d = pb_model(lat, pb.translational_symmetry())
    solver2d = pb_solver(mod2d)

    solver2d.set_wave_vector([kx, ky])
    Eks = solver2d.eigenvalues
    E1, E2 = Eks[band1], Eks[band2]
    vks = solver2d.eigenvectors
    v1 = vks[:, band1]
    solver2d.set_wave_vector([kx+qx, ky+qy])
    vkqs = solver2d.eigenvectors
    v2 = vkqs[:, band2]
    1e16*prefactor*(E1-E2)*sum(conj(v1).*v2)/sqrt(qx^2+qy^2)
end

function momentum(k::Vector{<:Real}, q::Vector{<:Real}, band1::Integer, lat::PyCall.PyObject)
    prefactor = (0.5*1e6/(3e18)^2)/(6.6e-16)
    kx, ky = k
    qx, qy = q
    mod2d = pb_model(lat, pb.translational_symmetry())
    solver2d = pb_solver(mod2d)

    solver2d.set_wave_vector([kx, ky])
    Eks = solver2d.eigenvalues
    E1 = Eks[band1]
    solver2d.set_wave_vector([kx+qx, ky+qy])
    Ekqs = solver2d.eigenvalues
    E2 = Ekqs[band1]
    1e16*prefactor*(E1-E2)/sqrt(qx^2+qy^2)

end