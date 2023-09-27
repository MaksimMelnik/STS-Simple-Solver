function [R_exch_data, Q] = ...
            R_exch(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, reaction)
% The second iteration of the universal function for exchange reactions
% M1 + M2 -> M3 + M4. 
% M1 and M3 are molecules; M2 and M4 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini).
% Electronic excitaion is not taken into account.
% 22.07.2023 Maksim Melnik

M3.num_vibr_levels=38;
M3.ev_0(1)=1.8843e-20;
M3.ev_i{1}=[0	3.72667168899985e-20	7.39742492804083e-20...	
    1.10122597171229e-19	1.45711760562462e-19	1.80741739454105e-19...
    2.15212533846160e-19	2.49124143738626e-19	2.82476569131503e-19...
    3.15269810024792e-19	3.47503866418492e-19	3.79178738312603e-19...
    4.10294425707125e-19	4.40850928602059e-19	4.70848246997404e-19...
    5.00286380893160e-19	5.29165330289328e-19	5.57485095185906e-19...
    5.85245675582896e-19	6.12447071480297e-19	6.39089282878110e-19...
    6.65172309776333e-19	6.90696152174968e-19	7.15660810074015e-19...
    7.40066283473472e-19	7.63912572373341e-19	7.87199676773621e-19...
    8.09927596674312e-19	8.32096332075414e-19	8.53705882976928e-19...
    8.74756249378853e-19	8.95247431281190e-19	9.15179428683937e-19...
    9.34552241587096e-19	9.53365869990666e-19	9.71620313894647e-19...
    9.89315573299040e-19	1.00645164820384e-18];
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
  kb = kf;
  kf_eq = reaction.A(T) * T^reaction.n(T)*exp(-reaction.E / k/T ); % m^3/sec
%   kf(reaction.index{1}, reaction.index{3}) = ...
%       kf(reaction.index{1}, reaction.index{3}) + reaction.A;
    kf = ...
      kf + kf_eq;
  dE = (repmat(M3.ev_i{1} + M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
- repmat((M1.ev_i{1} + M1.ev_0(1))', 1, M3.num_vibr_levels(1))) + ...
                                        (M1.diss_e(1) - M3.diss_e(1));
  if reaction.reverse
   Theta_r_M1 = M1.Be(1) * h * c / k;
   Z_rot_M1 = T ./ (M1.sigma .* Theta_r_M1);
   Theta_r_M3 = M3.Be(1) * h * c / k;
   Z_rot_M3 = T ./ (M3.sigma .* Theta_r_M3);
   Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) ...
        * (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 ...
                                                       * exp( dE / (k*T));
    % rate coefficient of backward (b) reaction
   kb = kf .* Kfb;
  end

  R_exch_data = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;
  dE_Q = dE + M3.ev_0(1) + M3.ev_i{1} - M1.ev_0(1) + M1.ev_i{1}';
  Q = sum(- R_exch_data .* dE_Q, 'all');
  %warning('recheck the energy in Q')
 case "ATn"
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)
 case "Arrhenius"
    coll.ArrA = reaction.A(T);
    coll.ArrN = reaction.n(T);
    coll.ArrE = reaction.E / k;   % in K
    [R_exch_data, Q] = R_exch_Arrhenius(M1, M2, M3, M4, ...
                                n_M1, n_M2, n_M3, n_M4, T, coll);
 case "Heaviside"
%   error("Exchange reactions of this type are still not " + ...
%                                         "implemented " + reaction.type)
    % constants
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
  coll.ArrA = reaction.A(T);
  coll.ArrN = reaction.n(T);
  coll.ArrE = reaction.E / k;   % in K
  kf = R_exch_Heaviside(M1, M2, M3, M4, ...
                                n_M1, n_M2, n_M3, n_M4, T, coll);
  dE = (repmat(M3.ev_i{1} + M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
    - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, M3.num_vibr_levels(1))) + ...
                                        (M1.diss_e(1) - M3.diss_e(1));

  Theta_r_M1 = M1.Be(1)*V_H*V_C/V_K;  Z_rot_M1 = T./(M1.sigma.*Theta_r_M1);
  Theta_r_M3 = M3.Be(1)*V_H*V_C/V_K;  Z_rot_M3 = T./(M3.sigma.*Theta_r_M3);

  Kfb = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 * ...
      exp( dE / (V_K*T));
  Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) * Kfb;
  %rate of backward (b) reaction
  kb = kf .* Kfb;
  R_exch_data = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;
  Q = sum(- R_exch_data .* dE, 'all');
 case "Heaviside, avg"
         % constants
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
  coll.ArrA = reaction.A(T);
  coll.ArrN = reaction.n(T);
  coll.ArrE = reaction.E / k;   % in K
  load("..\data\particles.mat", "NO");
  NO.num_elex_levels=1; 
  kf = R_exch_Heaviside(M1, M2, M3, M4, ...
                                n_M1, n_M2, n_M3, n_M4, T, coll);
  kf=sum(kf,2);
  dE = (repmat(M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
    - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, 1)) + ...
                                            (M1.diss_e(1) - M3.diss_e(1));

  Theta_r_M1 = M1.Be(1)*V_H*V_C/V_K;  Z_rot_M1 = T./(M1.sigma.*Theta_r_M1);
  Theta_r_M3 = M3.Be(1)*V_H*V_C/V_K;  Z_rot_M3 = T./(M3.sigma.*Theta_r_M3);
  
  Kfb = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 * ...
      exp( dE / (V_K*T));
  Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) * Kfb;
  %rate of backward (b) reaction
  kb = kf .* Kfb;
  R_exch_data = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;
  Q = sum(- R_exch_data .* dE, 'all');
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)    
end

end

function kf = ...
   R_exch_Heaviside(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, coll)
M3.num_vibr_levels=38;
M3.ev_0(1)=1.8843e-20;
M3.ev_i{1}=[0	3.72667168899985e-20	7.39742492804083e-20	1.10122597171229e-19	1.45711760562462e-19	1.80741739454105e-19	2.15212533846160e-19	2.49124143738626e-19	2.82476569131503e-19	3.15269810024792e-19	3.47503866418492e-19	3.79178738312603e-19	4.10294425707125e-19	4.40850928602059e-19	4.70848246997404e-19	5.00286380893160e-19	5.29165330289328e-19	5.57485095185906e-19	5.85245675582896e-19	6.12447071480297e-19	6.39089282878110e-19	6.65172309776333e-19	6.90696152174968e-19	7.15660810074015e-19	7.40066283473472e-19	7.63912572373341e-19	7.87199676773621e-19	8.09927596674312e-19	8.32096332075414e-19	8.53705882976928e-19	8.74756249378853e-19	8.95247431281190e-19	9.15179428683937e-19	9.34552241587096e-19	9.53365869990666e-19	9.71620313894647e-19	9.89315573299040e-19	1.00645164820384e-18];

    % constants
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

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

dE_M1 = coll.ArrE - ...
repmat((M1.ev_0(1) + M1.ev_i{1})'/V_K, 1, M3.num_vibr_levels(1)) + ...
repmat((M3.ev_0(1) + M3.ev_i{1})/V_K, M1.num_vibr_levels(1), 1);

EXP_M1 = exp(-dE_M1 .* heaviside(dE_M1) * (1/T));

kf_eq = coll.ArrA * T^coll.ArrN*exp(-coll.ArrE/T ); % m^3/sec

B_M1 = kf_eq / sum(EXP_M1 .* nM1i','all');

%rate of forward (f) reaction
kf = B_M1 * EXP_M1;

% if VDOP==0
% kf=sum(kf,2);
% end
% 
% if VDOP==0
% %ratio of rate of backward reaction to rate of forward reaction
% dE = (repmat(M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
%     - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, 1)) + ...
%                                             (M1.diss_e(1) - M3.diss_e(1));
% %dE=(M1.diss_e(1) - M3.diss_e(1));
% else
% dE = (repmat(M3.ev_i{1} + M3.ev_0(1), M1.num_vibr_levels(1), 1) ...
% - repmat((M1.ev_i{1} + M1.ev_0(1))', 1, M3.num_vibr_levels(1))) + ...
%                                         (M1.diss_e(1) - M3.diss_e(1));
% end
% 
% Kfb = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 * ...
%  exp( dE / (V_K*T));
% Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) * Kfb;
% 
% %rate of backward (b) reaction
% kb = kf .* Kfb;
% RExch1 = n_M3' * n_M4 .* kb  -  n_M1 * n_M2 .* kf;
% Q = sum(- RExch1 .* dE, 'all');
end



function [RExch1, Q] = ...
   R_exch_Arrhenius(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, coll)
    % constants
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

% parameters in Arrhenius law for reactuib in structure coll: ArrA, ArrN,
% ArrE
%(kd_eq=A*T^N*exp(-E/T))

Theta_r_M1 = M1.Be(1)*V_H*V_C/V_K;   Z_rot_M1 = T./(M1.sigma.*Theta_r_M1);
Theta_r_M3 = M3.Be(1)*V_H*V_C/V_K;   Z_rot_M3 = T./(M3.sigma.*Theta_r_M3);

kf_eq = coll.ArrA * T^coll.ArrN*exp(-coll.ArrE/T ); % m^3/sec

%rate of forward (f) reaction
kf = kf_eq;

%ratio of rate of backward reaction to rate of forward reaction
dE = (M1.diss_e(1) - M3.diss_e(1));
Kfb = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5 * Z_rot_M1/Z_rot_M3 * ...
 exp( dE / (V_K*T));
Kfb = M1.s_e(1) * M2.s_e(1) / (M3.s_e(1) * M4.s_e(1)) * Kfb;


%rate of backward (b) reaction
kb = kf .* Kfb;
RExch1 = sum(n_M3)' * n_M4 .* kb  -  sum(n_M1) * n_M2 .* kf;
Q = sum(- RExch1 .* dE, 'all');
end