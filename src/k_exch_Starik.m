function kd = k_exch_Starik(T, kd_eq, reaction, Ps_r, Ps_p)
% Exchange rate coefficient using the Starik model.
% kd is the vector of rate coefficients.
% M is the dissociating molecule; T is the gas temperature;
% kd_eq is the equlibrium rate coefficient;
% U is the non-equilibrium parameter.
% 29.09.2023 Maksim Melnik

k = 1.380649e-23;
n_eq = 1;
form_e_r = 0;
for ind_r = 1:length(Ps_r)
 n_eq_temp = density_f_exc(T, 1, Ps_r{ind_r});
 n_eq_temp = ...
        reshape(n_eq_temp, [ones(1, ind_r - 1), length(n_eq_temp), 1]);
 n_eq = n_eq .* n_eq_temp;
 form_e_r = form_e_r + Ps_r{ind_r}.form_e;
end
form_e_p = 0;
for ind_p = 1:length(Ps_p)
 form_e_p = form_e_p + Ps_p{ind_p}.form_e;
end
Ea = reaction.E;
ET = form_e_p - form_e_r;
% Ep = Epr(Ps_p);
    % product energy
Ep = 0;
for ind_p = 1:length(Ps_p)
 P1 = Ps_p{ind_p};
 if P1.fr_deg_c > 3
  e_ev = zeros(1, sum(P1.num_vibr_levels(1:P1.num_elex_levels)));
  for ind_e = 1:P1.num_elex_levels
   e_ev(1 + sum(P1.num_vibr_levels(1:ind_e-1)):...
                    sum(P1.num_vibr_levels(1:ind_e))) = ...
                        P1.e_E(ind_e) + P1.ev_0(ind_e) + P1.ev_i{ind_e};
  end
 else
  e_ev = P1.e_E(1:P1.num_elex_levels);
 end
 e_ev = reshape(e_ev, [ones(1, ind_p - 1), length(e_ev), 1]);
 Ep = Ep + e_ev;
end
% Er = Epr(Ps_r);
    % reagent energy
Er = 0;
for ind_r = 1:length(Ps_r)
 P1 = Ps_r{ind_r};
 if P1.fr_deg_c > 3
  e_ev = zeros(1, sum(P1.num_vibr_levels(1:P1.num_elex_levels)));
  for ind_e = 1:P1.num_elex_levels
   e_ev(1 + sum(P1.num_vibr_levels(1:ind_e-1)):...
                    sum(P1.num_vibr_levels(1:ind_e))) = ...
                        P1.e_E(ind_e) + P1.ev_0(ind_e) + P1.ev_i{ind_e};
  end
 else
  e_ev = P1.e_E(1:P1.num_elex_levels);
 end
 e_ev = reshape(e_ev, [ones(1, ind_p + ind_r - 1), length(e_ev), 1]);
 Er = Er + e_ev;
end
dE = max(Ea, ET + Ep) - Er;
expdE_Pneq = exp(- dE .* heaviside(dE) / k / T) .* n_eq;
sum_expdE_Pneq = sum(expdE_Pneq, 'all');
A = kd_eq * sum_dE_Pneq;
% kd = A * exp(- (dE * heaviside(dE)) / T);
% 
% k = 1.380649e-23;   % Boltzmann constant
% zero_template = zeros(1, sum(M.num_vibr_levels(1:M.num_elex_levels)));
% V  = zero_template; % non-equilibrium factor
% kd = zero_template; % diss coef
% for i = 1:M.num_elex_levels    % вот тут энергию e_i от нуля или от e_0?
%  if M.num_vibr_levels(i) > 1
%   %  ZvibT=sum(exp(-M.ev_i{i}/(k*T)));
%   %  ZvibU=sum(exp(M.ev_i{i}/(k*U(i))));
%   ZvibT = sum( exp(-(M.ev_0(i) + M.ev_i{i})/(k*T)) ); %vibrational statsum
%   ZvibU = sum( exp( (M.ev_0(i) + M.ev_i{i})/(k*U(i))) );
%   % V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=... old
%   %                exp((M.ev_i{i}+M.e_E(i))/k*(1/T+1/U(i)))*ZvibT/ZvibU;
%   V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...new
%         exp((M.ev_0(i) + M.ev_i{i})/k*(1/T + 1/U(i))) * ZvibT / ZvibU;
%   kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
%       V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) ...
%                                                                * kd_eq(i);
%  else
%   kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
%                                                                 kd_eq(i);
%  end
% end

end

% function Epr = Epr(Ps_r)
% % Energy calculation
% Epr = 0;
% for ind_r = 1:length(Ps_r)
%  P1 = Ps_r{ind_r};
%  if P1.fr_deg_c > 3
%   e_ev = zeros(1, sum(P1.num_vibr_levels(1:P1.num_elex_levels)));
%   for ind_e = 1:Ps_r{ind_r}.num_elex_levels
%    e_ev(1 + sum(P1.num_vibr_levels(1:ind_e-1)):...
%                     sum(P1.num_vibr_levels(1:ind_e))) = ...
%                         P1.e_E(ind_e) + P1.ev_0(ind_e) + P1.ev_i{ind_e};
%   end
%  else
%   e_ev = P1.e_E(1:P1.num_elex_levels);
%  end
%  e_ev = reshape(e_ev, [ones(1, ind_r - 1), length(e_ev), 1]);
%  Epr = Epr + e_ev;
% end
% end