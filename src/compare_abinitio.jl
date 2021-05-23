
"""
$(TYPEDSIGNATURES)
Compare a dft calculation from JDFTX with a tight binding calculation from PyBinding. 

"""
function compare_dft(dft_filebase::AbstractString, kpoints_file::AbstractString, models::Vector{<:PyCall.PyObject}, dft_dir::AbstractString="./")


end

"""
$(TYPEDSIGNATURES)
Compare a dft calculation from JDFTX with a tight binding calculation from PyBinding. 

"""
function compare_dft(dft_filebase::AbstractString, kpoints_file::AbstractString, model::PyCall.PyObject; 
    nbands::Integer=72, dft_dir::String="./", minband::Integer=1, maxband::Integer=1, kwargs...)
    Plots.plot()
    kpoints = Vector{Vector{Float64}}()
    for line in readlines(kpoints_file)[3:end]
        b1, b2 = model.lattice.reciprocal_vectors()
        push!(kpoints, sum(parse.(Float64, String.(split(line))[2:3]).*[b1, b2]) )
    end
    nkpts = length(kpoints)
    solver = pb_solver(model)
    bands = Vector{Vector{Float64}}()
    for kpoint in kpoints
        solver.set_wave_vector(kpoint)
        push!(bands, solver.eigenvalues)
    end
    bands_2d = zeros(length(bands), length(bands[1]) )
    for row in 1:length(bands)
        bands_2d[row, :]=bands[row]
    end
    np = pyimport("numpy")
    dftbands = np.reshape(np.fromfile(dft_dir*dft_filebase), (nkpts*2, nbands))*27.2
    a = Plots.plot!(bands_2d, color="pink", label="TB Model"; kwargs...)
    Plots.plot!(dftbands[1:nkpts, minband:maxband], label="DFT Spin Up", legend=false, color="blue"; kwargs...)
    Plots.plot!(dftbands[nkpts+1:2*nkpts, minband:maxband], label="DFT Spin Dn", color="red"; kwargs...)
    display(Plots.plot!(a, xticks=[],legend=false, size=(1000, 500)))
end

function compare_dft(dft_filebases::Vector{<:AbstractString}, kpoints_file::AbstractString, models::Vector{<:PyCall.PyObject}; 
    nbands::Integer=72, dft_dir::String="./", minband::Integer=1, maxband::Integer=1, kwargs...)

    Plots.plot()
    for (dft_filebase, model) in zip(dft_filebases, models)
        compare_dft(dft_filebase, kpoints_file, model; 
            nbands=nbands, dft_dir=dft_dir, minband=minband, maxband=maxband, kwargs...)
    end

end