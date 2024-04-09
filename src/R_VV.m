function [R, Q] = R_VV(M1, n_m, M2, n_M2, T, ind_e1, ind_e2, model)
% Universal function for relaxation parts R_VV calculation.
% R is the relaxation terms vector, Q is the flux.
% M1 is the first molecule, n_m is the array of M1's number
% densities (n_i), M2 is the collision partner, n_M2 is the number density
% of M2, T is the gas temperature, ind_e1 is the electronic level of M1,
% ind_e2 is the electronic level of M2, model (optional) is the VT model.
% 20.01.2023 Maksim Melnik

if nargin<8
    model='SSH';
end

k = 1.380649e-23;

switch model
    case 'FHO'
        k_down=kvv_fho_old(T, M1, M2, ind_e1, ind_e2);
    case 'FHO-FR'
        k_down=kvv_fho_old(T, M1, M2, ind_e1, ind_e2);
    case 'SSH'
        k_down=kvv_ssh(T, M1, M2, ind_e1, ind_e2, 1);
    case 'Guerra'
        k_down = kvv_ssh(T, M1, M2, ind_e1, ind_e2, 1) / 3e1;
end
dE1 = (M1.ev_i{ind_e1}(2:end) - M1.ev_i{ind_e1}(1:end-1))';
dE2 = M2.ev_i{ind_e2}(2:end) - M2.ev_i{ind_e2}(1:end-1);
k_up = k_down .* exp((dE2 - dE1)/k/T);      % detailed balance
R = zeros(M1.num_vibr_levels(ind_e1), M2.num_vibr_levels(ind_e2));
n_M2 = n_M2';
R_down = n_m(2:end)   .* n_M2(1:end-1) .* k_down;
R_up   = n_m(1:end-1) .* n_M2(2:end)   .* k_up;
R(2:end, 1:end-1) = R_up - R_down;
R(1:end-1, 2:end) = R(1:end-1, 2:end) - R_up + R_down;
Q = sum(R_down .* (dE1 - dE2) + R_up .* (dE2 - dE1), 'all');
end