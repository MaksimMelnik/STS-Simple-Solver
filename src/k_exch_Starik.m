function [kd, dE_fb] = k_exch_Starik(T, kd_eq, reaction, Ps_r, Ps_p)
% Exchange rate coefficient using the Starik model [1].
% kd is the vector of rate coefficients.
% M is the dissociating molecule; T is the gas temperature;
% kd_eq is the equlibrium rate coefficient;
% U is the non-equilibrium parameter.
% [1] B.I. Loukhovitski, A.M. Starik / Chemical Physics 360 (2009) 18–26
% 29.09.2023 Maksim Melnik

k = 1.380649e-23;
n_eq = 1;
form_e_r = 0;
form_e_p = 0;
    % reagents
Er = 0;
for ind_r = 1:length(Ps_r)
 P1 = Ps_r{ind_r};
 alpha_c = reaction.index{ind_r}{1};
 n_eq_temp = distribution_Boltzmann(T, 1, P1, alpha_c);
 n_eq_temp = ...
        reshape(n_eq_temp, [ones(1, ind_r - 1), length(n_eq_temp), 1]);
 n_eq = n_eq .* n_eq_temp;
 form_e_r = form_e_r + P1.form_e;
 if P1.fr_deg_c > 3
  e_ev = zeros(1, sum(P1.num_vibr_levels(alpha_c)));
  for ind_e = 1:length(alpha_c)
   e_ev(1 + sum(P1.num_vibr_levels(alpha_c(1:ind_e-1))): ...
                    sum(P1.num_vibr_levels(alpha_c(1:ind_e)))) = ...
                    P1.e_E(alpha_c(ind_e)) + ...
                        P1.ev_0(alpha_c(ind_e)) + P1.ev_i{alpha_c(ind_e)};
  end
 else
  e_ev = P1.e_E(alpha_c);
 end
 e_ev = reshape(e_ev, [ones(1, ind_r - 1), length(e_ev), 1]);
 Er = Er + e_ev;
end
    % products
Ep = 0;
for ind_p = 1:length(Ps_p)
 P1 = Ps_p{ind_p};
 form_e_p = form_e_p + P1.form_e;
 alpha_c = reaction.index{ind_r + ind_p}{1};
 if P1.fr_deg_c > 3
  e_ev = zeros(1, sum(P1.num_vibr_levels(alpha_c)));
  for ind_e = 1:length(alpha_c)
   e_ev(1 + sum(P1.num_vibr_levels(alpha_c(1:ind_e-1))): ...
                    sum(P1.num_vibr_levels(alpha_c(1:ind_e)))) = ...
                        P1.e_E(alpha_c(ind_e)) + ...
                        P1.ev_0(alpha_c(ind_e)) + P1.ev_i{alpha_c(ind_e)};
  end
 else
  e_ev = P1.e_E(alpha_c);
 end
 e_ev = reshape(e_ev, [ones(1, ind_r + ind_p - 1), length(e_ev), 1]);
 Ep = Ep + e_ev;
end
Ea = reaction.E;
ET = form_e_p - form_e_r;
% dE = (ET + Ep) - Er;
dE = max(Ea, ET + Ep) - Er;
exp_dE_Heav = exp(- dE .* heaviside(dE) / k / T);
expdE_Pneq = exp_dE_Heav .* n_eq;
sum_expdE_Pneq = sum(expdE_Pneq, 'all');
A = kd_eq / sum_expdE_Pneq;
kd = A * exp_dE_Heav;       % the coefficient
dE_fb = ET + Ep - Er;       % for the detailed balance
end