"""
$(TYPEDSIGNATURES)
Plots the defect atom (assumed to be in first position) with different color from the rest and also plots hoppings in such a way as to indicate differences in energies
"""
function plot_defectmodel(defect_model::PyCall.PyObject, natoms::Integer)
    defect_model.plot(hopping=Dict("width" =>6, "cmap"=>"auto"), site=Dict("cmap"=>[(i==1 ? "purple" : "red") for i in 1:natoms])); plt.show()
    plt.savefig("defectmodel.png")
    plt.close()
end