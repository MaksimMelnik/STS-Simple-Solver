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
% c = 299792458;            % speed of light
% h = 6.626070041e-34;      % Plank constant, J*sec

index = cell(1, 4);       % indexes from the reaction structure
Ms = {M1, M2, M3, M4};
i_sum = cell(1, 4);     % indexes from the 0th lvl of the ground el. state
E_t   = cell(1, 4);       % total energy
for ind = 1:4
 index{ind} = reaction.index{ind};
 if class(index{ind}{2}) == "string"
  index{ind}{2} = 1:Ms{ind}.num_vibr_levels(index{ind}{1});
 end
 i_sum{ind} = ...
        sum(Ms{ind}.num_vibr_levels(1:index{ind}{1}-1)) + index{ind}{2};
 E_t{ind}   = Ms{ind}.form_e + Ms{ind}.e_E(index{ind}{1});
 if Ms{ind}.fr_deg_c > 3
  E_t{ind}  = E_t{ind} + Ms{ind}.ev_0(index{ind}{1}) ...
                            + Ms{ind}.ev_i{index{ind}{1}}(index{ind}{2});
 end
end

kf = zeros(length(index{1}{2}), length(index{2}{2}), ...
                                length(index{3}{2}), length(index{4}{2}));
	% rate coefficient of backward (b) reaction
kb = kf;

switch reaction.type
 case "const"
  kf = kf + reaction.A;
  dE_fb = M3.form_e + M4.form_e - M1.form_e - M2.form_e;
 case "ATn"
  kf = kf + reaction.A * T ^ reaction.n;
  dE_fb = M3.form_e + M4.form_e - M1.form_e - M2.form_e;
 case "A(T/d_T)^n"
  kf = kf + reaction.A * (T / reaction.d_T) ^ reaction.n;
  dE_fb = M3.form_e + M4.form_e - M1.form_e - M2.form_e;
 case "Arrhenius"
  kf = kf + reaction.A * T ^ reaction.n * exp(- reaction.E / k / T);
  dE_fb = M3.form_e + M4.form_e - M1.form_e - M2.form_e;
 case "Heaviside"
     error("rewrite")
  [kf, dE_fb] = R_exch_Heaviside(M1, M3, T, reaction);
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)    
end
    % if the reaction proceeds in the opposite direction
if ~reaction.direction_forward
 kb = kf;
 kf = kf * 0;
end
if reaction.reverse     % if backward reaction included
    error('rewrite')
 Theta_r_M1 = M1.Be(1) * h * c / k;
 Z_rot_M1 = T ./ (M1.sigma .* Theta_r_M1);  % statistical rotational sum
 Theta_r_M3 = M3.Be(1) * h * c / k;
 Z_rot_M3 = T ./ (M3.sigma .* Theta_r_M3);
 Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) ...
        * (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 ...
                                                   * exp( dE_fb / (k*T));
 if reaction.direction_forward
  kb = kf .* Kfb;
 else
  kf = kb ./ Kfb;
 end
end
R_exch_data = ...
            zeros(length(n_M1), length(n_M2), length(n_M3), length(n_M4));
R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}) = ... 
    reshape(n_M3(i_sum{3}), 1, 1, []) .* ...
                        reshape(n_M4(i_sum{4}), 1, 1, 1, []) .* kb  ...
                                - n_M1(i_sum{1}) .* n_M2(i_sum{2})' .* kf;
    % energy change in the reaction (after - before)
dE_Q = reshape(E_t{3}, 1, 1, []) + reshape(E_t{4}, 1, 1, 1, []) ...
                                                    - E_t{1}' - E_t{2};
Q = sum(R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}) .* dE_Q, ...
                                                                'all');
end

function [kf, dE_fb] = R_exch_Heaviside(M1, M3, T, reaction)
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