function [R_exch_data, Q, R_exch_data_1, R_exch_data_2, ...
    R_exch_data_3, R_exch_data_4] = ...
                        R_exch(Ms, n_M1, n_M2, n_M3, n_M4, T, reaction)
% The third iteration of the universal function for exchange reactions
% M1 + M2 -> M3 + M4. 
% Ms are involved particles; M2 and M4 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini).
% Electronic excitaion is not taken into account.
% 22.07.2023 Maksim Melnik

    % constants
k = 1.380649e-23;         % Boltzmann constant, J/K
c = 299792458;            % speed of light
h = 6.626070041e-34;      % Plank constant, J*sec

num_particles = length(Ms);
index = cell(1, num_particles);     % indexes from the reaction structure
    % indexes from the 0th lvl of the ground el. state
i_sum = cell(1, num_particles); 
E_t   = cell(1, num_particles);     % total energy
for ind = 1:num_particles
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
    % initialization of kf with a proper size
kf = zeros(length(index{1}{2}), length(index{2}{2}), ...
                                length(index{3}{2}), length(index{4}{2}));
if ~reaction.direction_forward  % transposition for an inverted direction
 kf = permute(kf, [3 4 1 2]);     
end
	% rate coefficient of backward (b) reaction
kb = kf;


switch reaction.neq_model    % choosing reaction type
 case {"const", "ATn", "Arrhenius", "equal", "equilibrium"}
  kf = kf + k_equilibrium(reaction, T);
  dE_fb = Ms{3}.form_e + Ms{4}.form_e - Ms{1}.form_e - Ms{2}.form_e;
 case "Heaviside"
  coll.ArrA = reaction.A(T);
  coll.ArrN = reaction.n(T);
  coll.ArrE = reaction.E / k;   % in K
  Ms2 = Ms;
  n_M12 = n_M1;
  n_M22 = n_M2;
  n_M32 = n_M3;
  n_M42 = n_M4;
  if ~reaction.direction_forward
   Ms2 = {Ms{3}, Ms{4}, Ms{1}, Ms{2}};
   n_M12 = n_M3;
   n_M22 = n_M4;
   n_M32 = n_M1;
   n_M42 = n_M2;
  end
  kf = k_exch_Heaviside(Ms2{1}, Ms2{2}, Ms2{3}, Ms2{4}, ...
                                n_M12, n_M22, n_M32, n_M42, T, coll);
  % kf = k_exch_Heaviside(Ms{1}, Ms{3}, T, coll);
  size_kf = size(kf);
  kf = reshape(kf, size_kf(1), 1, size_kf(2));
  dE_fb = (repmat(Ms{3}.ev_i{1} + Ms{3}.ev_0(1), ...
                                        Ms{1}.num_vibr_levels(1), 1) ...
    - repmat((Ms{1}.ev_i{1} + Ms{1}.ev_0(1))', 1, ...
                        Ms{3}.num_vibr_levels(1))) + ...
                                    (Ms{1}.diss_e(1) - Ms{3}.diss_e(1));
  size_dE_fb = size(dE_fb);
  dE_fb = reshape(dE_fb, size_dE_fb(1), 1, size_dE_fb(2));
 case "Heaviside, avg"
  coll.ArrA = reaction.A(T);
  coll.ArrN = reaction.n(T);
  coll.ArrE = reaction.E / k;   % in K
  load("..\data\particles.mat", "NO");
  NO.num_elex_levels=1; 
  kf = k_exch_Heaviside(Ms{1}, Ms{2}, Ms{3}, Ms{4}, ...
                                n_M1, n_M2, n_M3, n_M4, T, coll);
  kf=sum(kf,2);
  dE_fb = (repmat(Ms{3}.ev_0(1), Ms{1}.num_vibr_levels(1), 1) ...
    - repmat((Ms{1}.ev_i{1} + Ms{1}.ev_0(1))', 1, 1)) + ...
                                    (Ms{1}.diss_e(1) - Ms{3}.diss_e(1));
 case {"Starik_test", "A(T/d_T)^n"}
  kd_eq = k_equilibrium(reaction, T);
  Ps_r = {Ms{1}, Ms{2}};
  Ps_p = {Ms{3}, Ms{4}};
  reaction2 = reaction;
  if ~reaction.direction_forward
   [Ps_p, Ps_r] = deal(Ps_r, Ps_p);     % swap
   reaction2.index{1} = reaction.index{3};
   reaction2.index{2} = reaction.index{4};
   reaction2.index{3} = reaction.index{1};
   reaction2.index{4} = reaction.index{2};
  end
  [kf, dE_fb] = k_exch_Starik(T, kd_eq, reaction2, Ps_r, Ps_p);
  % U = 3 * T;
  % [kf, dE_fb] = k_exch_Savelev(T, kd_eq, reaction2, Ps_r, Ps_p, U);
  if ~reaction.direction_forward
   dE_fb = - permute(dE_fb, [3 4 1 2]);
  end
 case "Starik_test_on_T"
  kd_eq = k_equilibrium(reaction, T);
  Ps_r = {Ms{1}, Ms{2}};
  Ps_p = {Ms{3}, Ms{4}};
  reaction2 = reaction;
  reaction2.E = reaction.E(T);
  if ~reaction.direction_forward
   [Ps_p, Ps_r] = deal(Ps_r, Ps_p);     % swap
   reaction2.index{1} = reaction.index{3};
   reaction2.index{2} = reaction.index{4};
   reaction2.index{3} = reaction.index{1};
   reaction2.index{4} = reaction.index{2};
  end
  [kf, dE_fb] = k_exch_Starik(T, kd_eq, reaction2, Ps_r, Ps_p);
  if ~reaction.direction_forward
   dE_fb = - permute(dE_fb, [3 4 1 2]);
  end
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)
end


    % if the reaction proceeds in the opposite direction
if ~reaction.direction_forward
 kb = permute(kf, [3 4 1 2]);     % transposition
 kf = kb * 0;
end
if reaction.reverse     % if backward reaction included
 Z_rot = zeros(1, length(Ms)) + 1;
 for ind = 1:length(Ms)
  if Ms{ind}.fr_deg_c > 3
   Z_rot(ind) = T / ... % statistical rotational sum
                (Ms{ind}.Be(index{ind}{1}) * h * c / k * Ms{ind}.sigma);
  end
 end
 Kfb = Ms{1}.s_e(index{1}{1}) * Ms{2}.s_e(index{2}{1}) / ...
     (Ms{3}.s_e(index{3}{1}) * Ms{4}.s_e(index{4}{1})) * ...
        (Ms{1}.mass * Ms{2}.mass / (Ms{3}.mass * Ms{4}.mass))^1.5 * ...
        prod(Z_rot(1:2)) / prod(Z_rot(3:end)) * exp( dE_fb / (k*T));
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
switch reaction.neq_model
 case {"equilibrium", "equal"}
  kf_eq = k_equilibrium(reaction, T);
  kb_eq = 0;
  if ~reaction.direction_forward
   kb_eq = kf_eq;
   kf_eq = 0;
  end
  if reaction.reverse     % if backward reaction included
   Z_rot = zeros(1, length(Ms)) + 1;
   for ind = 1:length(Ms)
    if Ms{ind}.fr_deg_c > 3
     Z_rot(ind) = T / ... % statistical rotational sum
                    (Ms{ind}.Be(index{ind}{1}) * h * c / k * Ms{ind}.sigma);
    end
   end
   Kfb_eq = Ms{1}.s_e(index{1}{1}) * Ms{2}.s_e(index{2}{1}) / ...
     (Ms{3}.s_e(index{3}{1}) * Ms{4}.s_e(index{4}{1})) * ...
        (Ms{1}.mass * Ms{2}.mass / (Ms{3}.mass * Ms{4}.mass))^1.5 * ...
            prod(Z_rot(1:2)) / prod(Z_rot(3:end)) * exp( dE_fb / (k*T));
   if reaction.direction_forward
    kb_eq = kf_eq .* Kfb_eq;
   else
    kf_eq = kb_eq ./ Kfb_eq;
   end
  end
  R_exch_data_1 = zeros(length(n_M1), 1);
  R_exch_data_2 = zeros(length(n_M2), 1);
  R_exch_data_3 = zeros(length(n_M3), 1);
  R_exch_data_4 = zeros(length(n_M4), 1);
  R_exch_data_1(i_sum{1}) = length(i_sum{2}) * sum(n_M3(i_sum{3}) ...
                                .* n_M4(i_sum{4})', 'all') * kb_eq ...
        - n_M1(i_sum{1}) * sum(n_M2(i_sum{2})) * length(i_sum{3}) ...
                                            * length(i_sum{4}) * kf_eq;
  R_exch_data_2(i_sum{2}) = length(i_sum{1}) * sum(n_M3(i_sum{3}) ...
                                .* n_M4(i_sum{4})', 'all') * kb_eq ...
        - sum(n_M1(i_sum{1})) * n_M2(i_sum{2}) * length(i_sum{3}) ...
                                            * length(i_sum{4}) * kf_eq;
  R_exch_data_3(i_sum{3}) = length(i_sum{1}) * length(i_sum{2}) ...
                    * n_M3(i_sum{3}) * sum(n_M4(i_sum{4})) * kb_eq ...
        - length(i_sum{4}) ...
                * sum(n_M1(i_sum{1}) .* n_M2(i_sum{2})', 'all') * kf_eq;
  R_exch_data_4(i_sum{4}) = length(i_sum{1}) * length(i_sum{2}) ...
                    * sum(n_M3(i_sum{3})) * n_M4(i_sum{4}) * kb_eq ...
        - length(i_sum{3}) ...
                * sum(n_M1(i_sum{1}) .* n_M2(i_sum{2})', 'all') * kf_eq;
  % Q2 = kb_eq * (length(i_sum{1}) * length(i_sum{2}) * ...
  %  sum(n_M3(i_sum{3}) .* n_M4(i_sum{4})' .* (E_t{3}' + E_t{4}), 'all') ...
  %  - (length(i_sum{2}) * sum(E_t{1}) + length(i_sum{1}) * sum(E_t{2})) ...
  %                   * sum(n_M3(i_sum{3}) .* n_M4(i_sum{4})', 'all')) ...
  %     - 0;
    otherwise
  R_exch_temp = R_exch_data;
  R_exch_temp4 = sum(R_exch_temp, 4);     % optimization, 
  R_exch_temp34 = sum(R_exch_temp4, 3);   %   this dumb way is faster
  R_exch_temp234 = sum(R_exch_temp34, 2);
  R_exch_temp134 = sum(R_exch_temp34, 1)';
  R_exch_temp24 = sum(R_exch_temp4, 2);
  R_exch_temp124 = reshape(sum(R_exch_temp24, 1), [], 1);
  R_exch_data_1 = R_exch_temp234;
  R_exch_data_2 = R_exch_temp134;
  R_exch_data_3 = R_exch_temp124;
  R_exch_data_4 = reshape(sum(R_exch_temp, [1 2 3]), [], 1);
end
    % energy change in the reaction (after - before)
dE_Q = reshape(E_t{3}, 1, 1, []) + reshape(E_t{4}, 1, 1, 1, []) ...
                                                    - E_t{1}' - E_t{2};
Q = sum(R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}) .* dE_Q, ...
                                                                'all');
end


% function kf = k_exch_Heaviside(M1, M3, T, coll)
function kf = k_exch_Heaviside(M1, M2, M3, M4, n1, n2, n3, n4, T, coll)
% parameters in Arrhenius law for reactuib in structure coll: ArrA, ArrN,
% ArrE
%(kd_eq=A*T^N*exp(-E/T))
k = 1.380649e-23; % Boltzmann constant, J/K

%vibrational statistical sums
exp_M1 = exp(-M1.ev_i{1}/(k*T));
Zv_M1 = sum(exp_M1);

% reduced Boltzmann distribution
nM1i = exp_M1 / Zv_M1;

dE_M1 = coll.ArrE - ...
repmat((M1.ev_0(1) + M1.ev_i{1})'/k, 1, M3.num_vibr_levels(1)) + ...
repmat((M3.ev_0(1) + M3.ev_i{1})/k, M1.num_vibr_levels(1), 1);

EXP_M1 = exp(-dE_M1 .* heaviside(dE_M1) * (1/T));

kf_eq = coll.ArrA * T^coll.ArrN*exp(-coll.ArrE/T ); % m^3/sec

B_M1 = kf_eq / sum(EXP_M1 .* nM1i','all');

%rate of forward (f) reaction
kf = B_M1 * EXP_M1;
end