function [RExch1, Q] = ...
            R_exch(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, coll, VDOP)
% Function for exchange reaction M1 + M2 -> M3 + M4.
% M1 and M3 are molecules; M2 and M4 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% coll is the collision variable; VDOP is a switcher of Vibrational 
% deactivation of product (1 - taking into account vibrational activation
% of product of reaction, 0 - without vibr. activation of the product 
% of the reaction).
% Electronic excitaion is not taken into account.
% 28.06.2023 Denis Kravchenko

    % constants
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

if (M2.fr_deg_c + M4.fr_deg_c ~= 6) 
    error("M2 and M4 must be atoms");
end
% parameters in Arrhenius law for reactuib in structure coll: ArrA, ArrN,
% ArrE
%(kd_eq=A*T^N*exp(-E/T))

Theta_r_M1 = M1.Be(1)*V_H*V_C/V_K;   Z_rot_M1 = T./(M1.sigma.*Theta_r_M1);
Theta_r_M3 = M3.Be(1)*V_H*V_C/V_K;   Z_rot_M3 = T./(M3.sigma.*Theta_r_M3);

%vibrational statistical sums
exp_M1 = exp(-M1.ev_i{1}/(V_K*T));
Zv_M1 = sum(exp_M1);

% reduced Boltzmann distribution
nM1i = exp_M1 / Zv_M1;


% difference of energy, K
dE_M1 = coll.ArrE - ...
    repmat((M1.ev_0(1) + M1.ev_i{1})'/V_K, 1, M3.num_vibr_levels(1)) + ...
    repmat((M3.ev_0(1) + VDOP*M3.ev_i{1})/V_K, M1.num_vibr_levels(1), 1);

EXP_M1 = exp(-dE_M1 .* heaviside(dE_M1) * (1/T));

kf_eq = coll.ArrA * T^coll.ArrN*exp(-coll.ArrE/T ); % m^3/sec

B_M1 = kf_eq / sum(EXP_M1 .* nM1i','all');

%rate of forward (f) reaction
kf = B_M1 * EXP_M1;

%ratio of rate of backward reaction to rate of forward reaction
dE = (repmat(M3.ev_i{1}*VDOP + M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
    - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, M3.num_vibr_levels(1))) + ...
                                            (M1.diss_e(1) - M3.diss_e(1));
Kfb = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 * ...
 exp( dE / (V_K*T));
Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) * Kfb;

%rate of backward (b) reaction
kb = kf .* Kfb;
RExch1 = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;

Q = sum(- RExch1 .* dE, 'all');
end