function [R, Q] = R_rec_wall(M, n_a, T, kinetics)
% Universal function for calculation of R_rec recombination terms on 
% the wall of the tube with radius R. The equations were taken from [1].
% 2 M + wall -> M2(X, 0)
% R is the vector of recombination terms R_rec_wall;
% Q is the energy flux of recombination.
% M is the molecule under consideration; n_a is the atom M's number 
% density; T is the gas temperature;
% kinetics is the big structure with all kinetics.
% 22.06.2023 Maksim Melnik
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006

k = 1.380649e-23;
eVtoJ = 1.602176565e-19;            % eV to J
switch M.name
    case 'O'
        gamma = 2.5e-3;             % DN
        dE = 5.12 * eVtoJ;          % eV
    case 'N'
        gamma = 1e-4;               % DN
        dE = 9.76 * eVtoJ;          % eV
    otherwise
        error("There're no data for wall recombination for this atom")
end
v_th = sqrt( 2*k*T / M.mass);       % thermal average velocity of M
nu_w = gamma * v_th / 2 / kinetics.tube_R;
R = - n_a .* nu_w;
beta = 0.5;
Q = - sum(R) * dE * (1-beta);
end