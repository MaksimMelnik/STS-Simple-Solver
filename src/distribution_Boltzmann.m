function n_ai = distribution_Boltzmann(T, n1, P, alpha)
% Boltzmann distribution function.
% n_ai is the vector of number densities for electronic states a and
% vibrational states i.
% T is the gas temperature; 
% n1 is the number density of particle P (use n1 = 1 for dimentionless);
% P is the particle structure; 
% alpha is the vector of desired electron states into accounting.
% Maksim Melnik
% 03.10.2023

k = 1.380649e-23;
Zel = P.s_e(alpha) * exp(- P.e_E(alpha) / k / T)';
s_e = P.s_e(alpha);
if length(alpha) == 1
 Zel = 1;
 s_e = 1;
 P.e_E(alpha) = 0;
end
s_vibr = 1;
if P.fr_deg_c > 3
 n_ai = zeros(1, sum(P.num_vibr_levels(alpha)));
 for ind_e = 1:length(alpha)
  Zvibr = sum(exp(-(P.ev_i{alpha(ind_e)} + P.ev_0(alpha(ind_e))) / k / T));
  e_ai = P.e_E(alpha(ind_e)) + P.ev_0(alpha(ind_e)) + P.ev_i{alpha(ind_e)};
  n_ai(1+sum(P.num_vibr_levels(alpha(1:ind_e-1))):...
            sum(P.num_vibr_levels(alpha(1:ind_e)))) = ...
                s_e(ind_e) * s_vibr / Zvibr * exp(- e_ai / k /T);
 end
else
 n_ai = s_e .* exp(- P.e_E(alpha) / k / T);
end
n_ai = n_ai * n1 / Zel;
end