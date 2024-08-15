classdef const
  properties(Constant)
    a = 2.84; %lattice parameter for LiF
    b = 5.48/3; %projected distance between Mo and S atoms
    c = 3.16; %sulphur-sulphur distance
    d = 5.48; %length of Min's unit cell

    D = 7.63; %Data for LiF
    alpha = 1.1;
    z0 = 1.0;
    beta = 0.1;

    MoS2Depth = 13.6%from https://pubs.acs.org/doi/epdf/10.1021/acs.jpcc.8b12029
    zOffset = 1;

    rotatMat = [sqrt(3)/2 -1/2;1/2 sqrt(3)/2];
    phong = 7;
  end
end