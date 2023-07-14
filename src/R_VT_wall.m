function [R, Q] = R_VT_wall(M, n_m, T, kinetics)
% Universal function for calculation of R_VT deactivation terms on 
% the wall of the tube with radius R. The equations were taken from [1].
% R is the vector of relaxation terms R_VT_wall, Q is the VT energy flux.
% M is the molecule under consideration; n_m is the array of M's number
% densities (n_i); T is the gas temperature;
% kinetics is the big structure with all kinetics.
% 21.06.2023 Maksim Melnik
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006

k = 1.380649e-23;

R = zeros(sum(M.num_vibr_levels(1:M.num_elex_levels)), 1);
gamma_v = 4.5e-4;                   % DN
v_th = sqrt( 2*k*T / M.mass);       % thermal average velocity of M
nu_w = gamma_v * v_th / 2 / kinetics.tube_R;
R_down = n_m .* nu_w;
for ind_e=1:M.num_elex_levels
 R_down(sum(M.num_vibr_levels(1:(ind_e-1)))+1) = 0; % the 0th vibr level
                                                %  of each el. exc. state
end
R(1:end-1) = R_down(2:end);
R = R - R_down;
e_i = [M.ev_i{1:M.num_elex_levels}];
beta = 0.5;
Q = sum((e_i(2:end) - e_i(1:end-1))' .* R_down(2:end)) * (1-beta);
end