"""
$(TYPEDSIGNATURES)

`subsampling` : Determines how much of the Brillouin zone is sampled (relevant for cases in which only certain regions of the Brillouin zone give non-negligible values for the polarization)
"""
function impol_2d(qx::Real, qy::Real, μ::Real, lat::PyCall.PyObject; spin::Integer=2, histogram_width::Real = 10, 
    mesh::Int=10, offset::Vector{<:Real}=[0, 0], subsampling::Real=1)

    qx *= 10 
    qy *= 10 ##The user is expected to give wavevectors in inverse angstrom so this converts to inverse nanometers
    im_pols = zeros(histogram_width*30)
    b1, b2 = lat.reciprocal_vectors()
    bzone_area = abs(cross(b1, b2)[3])/100 ## Brillouin zone area in inverse angstrom squared
    mod2d = pb_model(lat, pb.translational_symmetry())
    solver2d = pb_solver(mod2d)
    solver2d.set_wave_vector([0, 0])
    GammaBands = solver2d.eigenvalues
    num_bands = length(GammaBands)
    println(num_bands)
    for (xiter, yiter) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ##We offset the sampling of the brillouin zone to relevant k vectors if necessary and account for subsampling
        kx, ky = xiter/(subsampling*mesh)*b1[1] + yiter/(subsampling*mesh)*b2[1]+offset[1], xiter/(subsampling*mesh)*b1[2] + yiter/(subsampling*mesh)*b2[2]+offset[2]
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
                ω>0 && (im_pols[round(Int, ω*histogram_width+1)] += spin/(2π)^2*π*histogram_width*(f2-f1)*overlap*(1/mesh)^2*bzone_area*(1/subsampling)^2)
            end
        end
    end
    return im_pols
end

"""
$(TYPEDSIGNATURES)
"""
function realeps_2d(qx::Real, qy::Real, ωs::Vector{<:Real}, μ::Real, lat::PyCall.PyObject; histogram_width::Real=10, mesh::Int=100,
    max_energy_integration::Real=3, offset::Vector{<:Real}=zeros(2), subsampling::Real=1, spin::Integer=2, kwargs...)

    q=abs(sqrt(qx^2+qy^2))
    im_pol = impol_2d(qx, qy, μ, lat, mesh=mesh, histogram_width=histogram_width, offset=offset, subsampling=subsampling, spin=spin)
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(30-1/histogram_width), im_pol)
    epsilons = zeros(length(ωs))
    ErrorAbs=1e-10
    for (index, ω) in enumerate(ωs)
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
        epsilons[index] = 1-90.5/q*pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
    end
    return epsilons 
end

"""
$(TYPEDSIGNATURES)
Pass a vector of tuples corresponding to the wavevectors as well as a vector corresponding to the energies at which epsilon is evaluated and get the 
nonlocal dielectric function as a result. 
"""
function realepses_2d(qs::Vector{<:Tuple{<:Real, <:Real}}, ωs::Vector{<:Real}, μ::Real, lat::PyCall.PyObject; histogram_width::Real=10, mesh::Int=100, max_energy_integration::Real=3, offset::Vector{<:Real}=zeros(2), subsampling::Real=1, spin::Integer=2, kwargs...)
    numqs = length(qs)
    numomegas = length(ωs)
    epsilonarray = zeros(numqs, numomegas)
    println("The number of qs at which epsilon is being evaluated: ", numqs)
    println("The number of ωs at which epsilon is being evaluated: ", numomegas)
    for (index, qtuple) in enumerate(qs)
        println("Index: ", index)
        qx, qy = qtuple
        epsilonarray[index, :] = realeps_2d(qx, qy, ωs, μ, lat, histogram_width = histogram_width, mesh = mesh, max_energy_integration = max_energy_integration, offset = offset, subsampling= subsampling, spin = spin; kwargs...)
    end
    return epsilonarray
end

