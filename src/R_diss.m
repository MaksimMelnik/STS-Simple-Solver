function [RD, Q] = R_diss(M, n_m, n_a1, n_a2, n_p, Coll, T, n0, ind_U, ...
                                                                ind_diss)
% The universal function for calculating dissociation relaxation terms 
% R_Diss according to Arrhenius, Marrone-Treanor and Aliat models 
% for electronically excited molecules.
% RD is the relaxation terms R_Diss, Q is the energy flux.
% M is the dissociating molecule, n_m is the vector of concentrations,
% n_a1 is the concentration of the first atom of M1, n_a2 is the second 
% atom of M1, n_p is the concentration of collision partner particle,
% Coll is the collision parameters variable, T is the gas temperature,
% n0 is the characteristic number density, 
% ind_U is the non-equilibrium U parameter switcher, 
% ind_Aliat is the dissociation model switcher.
% 22.10.2020 Maksim Melnik

    % constants
V_K = 1.380649e-23;     V_C = 299792458;    V_H = 6.626070041e-34;
if nargin < 9       % default diss model and parameter
   ind_U = '3T';
   ind_diss = 'MT';
end

% Equilibrium dissociation rate coefficient using Arrhenius law.
% Coll contains the information about the collision. ArrA and ArrN are
% parameters in the Arrhenius law.
% kd_eq = Coll.ArrA * T^Coll.ArrN *... old
%                             exp(-(M.diss_e + M.e_E)/(V_K*T)); % m^3/sec
kd_eq = Coll.ArrA * T^Coll.ArrN *... new
                            exp(-(M.diss_e)/(V_K*T)); % m^3/sec

    % non-equilibrium parameter
U = M.diss_e * 0;
switch ind_U
    case 'D/6k'
        U = (M.diss_e + M.e_E) / (V_K*6);
    case '3T'
        U = U + 3*T;
    case 'inf'
        U = U + Inf;
    case 'Savelev16'    % only for O2
        U = U + M.diss_e/V_K * (0.15 + T/20000);
end
switch ind_diss     % choosing dissociation model
 case "Aliat"
  kd = k_diss_Aliat(M, T, kd_eq, U);
 case "MT"          % Marrone-Treanor model
  kd = k_diss_Marrone_Treanor(M, T, kd_eq, U);
 case "Arrhenius"
  error("Arrhenius model is still not implemented.")
    case "Savelev"     % perspective Savelev's model
  error("Savelev's model is still not implemented.")
end

Theta_r = M.Be * V_H * V_C / V_K;
Z_rot = T ./ (M.sigma .* Theta_r);  % rotational statistical sum
zero_template = zeros(1, sum(M.num_vibr_levels(1:M.num_elex_levels)));
dE_OA = zero_template;  % dE for energy flux E
Kdr   = zero_template;  % detailed balance coefficient
for i = 1:M.num_elex_levels
 dE = M.e_E(i) + M.ev_0(i) + M.ev_i{i} + M.form_e - M.form_e_atoms_sum;
 Kdr(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
                Z_rot(i) * M.s_e(i)/M.mltpl_atoms_s_e * exp(-(dE)/V_K/T);
 dE_OA(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = dE;
end
Kdr = Kdr * ...
        (M.mass/M.mltpl_atoms_mass)^(3/2) * V_H^3 * (2*pi*V_K*T)^(-3/2);
kr = kd .* Kdr * n0;    % recombination coefficient
RD = (n_p * (n_a1 * n_a2 * kr - (n_m') .* kd))';    % R_diss
Q = dE_OA * RD;         % energy flux Q
end