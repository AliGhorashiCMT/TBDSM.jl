function graphene_bands()
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    graphene_mod = model(pb_graphene.monolayer(), pb.translational_symmetry())
    solver(graphene_mod).calc_bands(K1, Gamma, M, K2).plot()
end

function bilayer_graphene_bands()
    a_cc=pb_graphene.a_cc
    Gamma = [0, 0]
    K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]
    M = [0, 2*pi / (3*a_cc)]
    K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
    graphene_mod = model(pb_graphene.bilayer(), pb.translational_symmetry())
    solver(graphene_mod).calc_bands(K1, Gamma, M, K2).plot()
end