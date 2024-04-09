function [R, Q] = R_VT(M1, n_m, M2, n_M2, T, ind_e, model, coefs_for_polys)
% Universal function for calculation of R_VT relaxation terms.
% R is the vector of relaxation terms R_VT, Q is the VT energy flow.
% M1 is the first molecule, n_m is the array of M1's number
% densities (n_i), M2 is the collision partner, n_M2 is the number density
% of M2, T is the gas temperature, ind_e is the electronic level of M1,
% model (optional) is the VT model.
% 11.10.2020 Maksim Melnik

if nargin < 7
    model = 'FHO';
end

k = 1.380649e-23;

switch model
    case 'FHO-FR'
        k_down = kvt_fho_fr(T, M1, M2, coefs_for_polys);
    case 'FHO'
        k_down = kvt_fho_old(T, M1, M2, ind_e);
    case 'SSH'
        k_down = kvt_ssh(T, M1, M2, ind_e, 1);
    case 'Guerra'
        if M1.name=="N2" && M2.name=="O"
            k_down = kvt_Gordietst(T, M1, ind_e);
        else
            k_down = kvt_ssh(T, M1, M2, ind_e, 1);
        end
end
k_down = k_down';
k_up  = k_down.*exp((M1.ev_i{ind_e}(1:end-1)-M1.ev_i{ind_e}(2:end))/k/T)';
core  = 1:M1.num_vibr_levels(ind_e) - 1;
R_down = n_M2 * n_m(core+1) .* k_down;
R_up   = n_M2 * n_m(core)   .* k_up;
R = zeros(M1.num_vibr_levels(ind_e), 1);
R(core)   = R_down - R_up;
R(core+1) = R(core+1) - R_down + R_up;
Q = (M1.ev_i{ind_e}(core+1)-M1.ev_i{ind_e}(core)) * (R_down - R_up);
end