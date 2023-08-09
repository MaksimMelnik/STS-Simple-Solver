function [R_exch_data, Q] = R_exch_23(...
            M1, M2, M3, M4, M5, n_M1, n_M2, n_M3, n_M4, n_M5, T, reaction)
% The temporal universal function for exchange reactions
% M1 + M2 -> M3 + M4 + M5. 
% M1, M2 and M3 are molecules; M4 and M5 are currently atoms; 
% n_Mi are number densities of Mi; T is the gas temperature; 
% reaction is the reaction structure (see data/reactions_data_ini).
% 08.08.2023 Maksim Melnik

if ~(M1.fr_deg_c>3 && M2.fr_deg_c>3 && M3.fr_deg_c>3 && M4.fr_deg_c==3 ...
                                                        && M5.fr_deg_c==3)
    error("Some of the particles are unsuitable.");
end
i1 = reaction.index{1};
i2 = reaction.index{2};
i3 = reaction.index{3};
i4 = reaction.index{4};
i5 = reaction.index{5};
kf = zeros(length(n_M1), length(n_M2), length(n_M3));
% kf = zeros(M1.num_vibr_levels(i1{1}), M2.num_vibr_levels(i2{1}), ...
%                                             M3.num_vibr_levels(i3{1}));
i1sum = sum(M1.num_vibr_levels(1:i1{1}-1)) + i1{2};
i2sum = sum(M2.num_vibr_levels(1:i2{1}-1)) + i2{2};
i3sum = sum(M3.num_vibr_levels(1:i3{1}-1)) + i3{2};
switch reaction.type
 case "A(T/d_T)^n"
     % it has a potential for optimizaation
%   kf(i1{2}, i2{2}, i3{2}) = reaction.A * (T / reaction.d_T) ^ reaction.n;
  kf(i1sum, i2sum, i3sum) = reaction.A * (T / reaction.d_T) ^ reaction.n;
 otherwise
        error("Exchange reactions of this type are still not " + ...
            "implemented " + reaction.type)    
end
   
	% rate coefficient of backward (b) reaction
kb = kf * 0;
    % if the reaction prceeds in the opposite direction
if ~reaction.direction_forward
  error(["Backward reactions in M + M -> M + A + A exchange are still" ...
        " are still not implemented"])
end
if reaction.reverse     % if backward reaction included
   error(["Reverse reactions in M + M -> M + A + A exchange are still" ...
        " are still not implemented"])
end
R_exch_data = reshape(n_M3, 1, 1, []) * n_M4(i4) * n_M5(i5) .* kb  -  ...
                                                    n_M1 .* n_M2' .* kf;
    % energy change in the reaction (after - before)
dE_Q = M3.form_e + M3.e_E(i3{1}) + M3.ev_0(i3{1}) ...
                                + reshape(M3.ev_i{i3{1}}, 1, 1, []) ...
	   + M4.form_e + M4.e_E(i4) + M5.form_e + M5.e_E(i5) ...
       - M1.form_e - M1.e_E(i1{1}) - M1.ev_0(i1{1}) - M1.ev_i{i1{1}}' ...
       - M2.form_e - M2.e_E(i2{1}) - M2.ev_0(i2{1}) - M2.ev_i{i2{1}};
Q = sum(R_exch_data(i1sum, i2sum, i3sum) .* dE_Q(i1{2}, i2{2}, i3{2}), ...
                                                                'all');
end