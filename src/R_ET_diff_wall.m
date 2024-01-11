function [R_ET, Q] = R_ET_diff_wall(M, n, T, kinetics)
% Universal function for calculation of R_ET deactivation terms on 
% the wall of the tube with radius R. The equations were taken from [1].
% R_ET is the vector of relaxation terms R_ET_wall, Q is the ET energy flux
% M is the molecule under consideration; n is the array of all number
% densities (n_ci); T is the gas temperature;
% kinetics is the big structure with all kinetics.
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006
% 03.08.2023 Maksim Melnik

beta = 0.5;
M_e_i = cell(1, M.num_elex_levels);
diff_cs = kinetics.reactions('ET');
diff_c = diff_cs{kinetics.IndexOfMolecules(M.name)};
for ind_e = 1:M.num_elex_levels
 M_e_i{ind_e} = ...
     1+sum(M.num_vibr_levels(1:ind_e-1)):sum(M.num_vibr_levels(1:ind_e));
end
n_m  = n(kinetics.index{kinetics.IndexOfMolecules(M.name)});
R_ET = n_m * 0;
dE   = n_m * 0;
N_g  = sum(n);
for ind_e = 2:M.num_elex_levels
 DY_N = diff_c(ind_e).A * (T / diff_c(ind_e).d_T) ^ diff_c(ind_e).n / N_g; %/1e4;
 R_ET(M_e_i{ind_e}) = DY_N;
 dE(M_e_i{ind_e})   = M.e_E(ind_e) + M.ev_0(ind_e) + M.ev_i{ind_e};
end
Lambda  = kinetics.tube_R / 2.4;
R_ET    = - R_ET .* n_m / Lambda^2;
R_ET(1) = - sum(R_ET(2:end));

Q = sum(- R_ET .* dE * (1 - beta));
end