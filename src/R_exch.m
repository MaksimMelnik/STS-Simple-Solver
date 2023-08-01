function [R_exch_data, Q] = ...
            R_exch(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, reaction)
% The second iteration of the universal function for exchange reactions
% M1 + M2 -> M3 + M4. 
% M1 and M3 are molecules; M2 and M4 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini).
% Electronic excitaion is not taken into account.
% 22.07.2023 Maksim Melnik

    % constants
k = 1.380649e-23;         % Boltzmann constant, J/K
c = 299792458;            % speed of light
h = 6.626070041e-34;      % Plank constant, J*sec

if (M2.fr_deg_c + M4.fr_deg_c ~= 6) 
    error("M2 and M4 must be atoms.");
end

switch reaction.type
 case "const"
  kf = zeros(M1.num_vibr_levels(1), M3.num_vibr_levels(1));
  kf(reaction.index{1}, reaction.index{3}) = ...
                    kf(reaction.index{1}, reaction.index{3}) + reaction.A;
  dE_fb = M1.diss_e(1) - M3.diss_e(1);
 case "ATn"
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)
 case "Arrhenius"
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)
 case "Heaviside"
  [kf, dE_fb] = R_exch_Heaviside(M1, M2, M3, M4, T, reaction);
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)    
end

	% rate coefficient of backward (b) reaction
kb = zeros(M1.num_vibr_levels(1), M3.num_vibr_levels(1));
if reaction.reverse     % if backward reaction included
 Theta_r_M1 = M1.Be(1) * h * c / k;
 Z_rot_M1 = T ./ (M1.sigma .* Theta_r_M1);  % statistical rotational sum
 Theta_r_M3 = M3.Be(1) * h * c / k;
 Z_rot_M3 = T ./ (M3.sigma .* Theta_r_M3);
 Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) ...
        * (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 ...
                                                   * exp( dE_fb / (k*T));
 kb = kf .* Kfb;
end
R_exch_data = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;
    % energy change in the reaction (after - before)
dE_Q = M1.diss_e(1) - M3.diss_e(1) ...
                    + M3.ev_0(1) + M3.ev_i{1} - M1.ev_0(1) - M1.ev_i{1}';
Q = sum(R_exch_data .* dE_Q, 'all');
end

function [kf, dE_fb] = R_exch_Heaviside(M1, M2, M3, M4, T, reaction)
    % constants
V_K = 1.380649e-23;

VDOP = 1;
if length(reaction.index{3}) == 1
  VDOP = 0;
end
coll.ArrA = reaction.A(T);
coll.ArrN = reaction.n(T);
coll.ArrE = reaction.E / V_K;   % in K

% parameters in Arrhenius law for reactuib in structure coll: ArrA, ArrN,
% ArrE
%(kd_eq=A*T^N*exp(-E/T))

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
dE_fb = (repmat(M3.ev_i{1}*VDOP + M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
    - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, M3.num_vibr_levels(1))) + ...
                                            (M1.diss_e(1) - M3.diss_e(1));
end