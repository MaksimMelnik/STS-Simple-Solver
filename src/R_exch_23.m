function [R_exch_data, Q] = R_exch_23(...
            M1, M2, M3, M4, M5, n_M1, n_M2, n_M3, n_M4, n_M5, T, reaction)
% The universal function for exchange reactions M1 + M2 -> M3 + M4 + M5. 
% M1, M2 and M3 are molecules; M4 and M5 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini).
% 08.08.2023 Maksim Melnik

num_particles = 5;
index = cell(1, num_particles);      % indexes from the reaction structure
Ms = {M1, M2, M3, M4, M5};
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
kb2 = kf;

switch reaction.type
 case "const"
  kf = kf + reaction.A;
 case "A(T/d_T)^n"
  kf = kf + reaction.A * (T / reaction.d_T) ^ reaction.n;
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
   error(["Reverse reactions in M + M -> M + A + A exchange are still" ...
        " are still not implemented"])
end
R_exch_data = zeros(length(n_M1), length(n_M2), length(n_M3), ...
                                            length(n_M4), length(n_M5));
R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}, i_sum{5}) = ... 
    reshape(n_M3(i_sum{3}), 1, 1, []) .* ...
                    reshape(n_M4(i_sum{4}), 1, 1, 1, []) .* ...
                    reshape(n_M5(i_sum{5}), 1, 1, 1, 1, []) .* kb2  ...
                            - n_M1(i_sum{1}) .* n_M2(i_sum{2})' .* kf;
    % energy change in the reaction (after - before)
dE_Q = reshape(E_t{3}, 1, 1, []) + reshape(E_t{4}, 1, 1, 1, []) ...
                    + reshape(E_t{5}, 1, 1, 1, 1, []) - E_t{1}' - E_t{2};
Q = sum(R_exch_data(i_sum{1}, i_sum{2}, i_sum{3}, i_sum{4}, i_sum{5}) ...
                                                        .* dE_Q, 'all');
end