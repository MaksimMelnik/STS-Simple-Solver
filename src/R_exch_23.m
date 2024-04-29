function [R_exch_data, Q] = ...
            R_exch_23(Ms, n_M1, n_M2, n_M3, n_M4, n_M5, T, reaction, n0)
% The universal function for exchange reactions M1 + M2 -> M3 + M4 + M5. 
% Ms are involved particles; M4 and M5 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini);
% n0 is the characteristic number density.
% 08.08.2023 Maksim Melnik

    % constants
k = 1.380649e-23;         % Boltzmann constant, J/K
c = 299792458;            % speed of light
h = 6.626070041e-34;      % Plank constant, J*sec

num_particles = length(Ms);
index = cell(1, num_particles);      % indexes from the reaction structure
% indexes from the 0th lvl of the ground el. state
i_sum = cell(1, num_particles);     
E_t   = cell(1, num_particles);       % total energy
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

kf = zeros(length(index{1}{2}), length(index{2}{2}), ...
        length(index{3}{2}), length(index{4}{2}), length(index{5}{2}));
	% rate coefficient of backward (b) reaction
kb = kf;

switch reaction.type
 case "const"
  kf = kf + reaction.A;
  dE_fb = Ms{3}.form_e + Ms{4}.form_e + Ms{5}.form_e ...
                                            - Ms{1}.form_e - Ms{2}.form_e;
 % case "A(T/d_T)^n"
 %  kf = kf + reaction.A * (T / reaction.d_T) ^ reaction.n;
 %  dE_fb = Ms{3}.form_e + Ms{4}.form_e + Ms{5}.form_e ...
 %                                            - Ms{1}.form_e - Ms{2}.form_e;
 case {"Starik_test", "A(T/d_T)^n"}
  kd_eq = reaction.A * (T / reaction.d_T) ^ reaction.n * ...
                                                exp(- reaction.E / k / T);
  Ps_r = {Ms{1}, Ms{2}};
  Ps_p = {Ms{3}, Ms{4}, Ms{5}};
  reaction2 = reaction;
  if ~reaction.direction_forward
      error("isn't tested yet")
   [Ps_p, Ps_r] = deal(Ps_r, Ps_p);     % swap
   reaction2 = reaction;
   reaction2.index{1} = reaction.index{3};
   reaction2.index{2} = reaction.index{4};
   reaction2.index{3} = reaction.index{1};
   reaction2.index{4} = reaction.index{2};
  end
  [kf, dE_fb] = k_exch_Starik(T, kd_eq, reaction2, Ps_r, Ps_p);
  U = 3 * T;
  % [kf, dE_fb] = k_exch_Savelev(T, kd_eq, reaction2, Ps_r, Ps_p, U);
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented: " + reaction.type)    
end
   
    % if the reaction prceeds in the opposite direction
if ~reaction.direction_forward
  error(["Backward reactions in M + M -> M + A + A exchange are still" ...
        " are still not implemented"])
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
     (Ms{3}.s_e(index{3}{1}) * Ms{4}.s_e(index{4}{1}) * ...
                                            Ms{5}.s_e(index{5}{1})) * ...
        (Ms{1}.mass * Ms{2}.mass / ...
                        (Ms{3}.mass * Ms{4}.mass * Ms{5}.mass))^1.5 * ...
        prod(Z_rot(1:2)) / prod(Z_rot(3:end)) * exp( dE_fb / (k*T)) * ...
        h ^ 3 * (2 * pi * k * T) ^ (-1.5);
 if reaction.direction_forward
  kb = kf .* Kfb * n0;
 else
  kf = kb ./ Kfb / n0;
 end
   % error(["Reverse reactions in M + M -> M + A + A exchange are still" ...
   %      " are still not implemented"])
end
R_exch_data = zeros(length(n_M1), length(n_M2), length(n_M3), ...
                                            length(n_M4), length(n_M5));
R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}, i_sum{5}) = ... 
    reshape(n_M3(i_sum{3}), 1, 1, []) .* ...
                    reshape(n_M4(i_sum{4}), 1, 1, 1, []) .* ...
                    reshape(n_M5(i_sum{5}), 1, 1, 1, 1, []) .* kb  ...
                            - n_M1(i_sum{1}) .* n_M2(i_sum{2})' .* kf;
    % energy change in the reaction (after - before)
dE_Q = reshape(E_t{3}, 1, 1, []) + reshape(E_t{4}, 1, 1, 1, []) ...
                    + reshape(E_t{5}, 1, 1, 1, 1, []) - E_t{1}' - E_t{2};
Q = sum(R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}, i_sum{5}) ...
                                                        .* dE_Q, 'all');
end