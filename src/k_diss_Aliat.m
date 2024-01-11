function kd = k_diss_Aliat(M, T, kd_eq, U)
% Dissociation coefficient using the Aliat model.
% kd is the vector of rate coefficients.
% M is the dissociating molecule; T is the gas temperature;
% kd_eq is the equlibrium rate coefficient;
% U is the non-equilibrium parameter.
% 24.10.2020 Maksim Melnik

V_K = 1.380649e-23;     % Boltzmann constant
Zel = M.s_e(1:M.num_elex_levels) * ...  electronic stat sum
                            exp( -M.e_E(1:M.num_elex_levels)' / (V_K*T));
divsum = 0;
zero_template = zeros(1, sum(M.num_vibr_levels(1:M.num_elex_levels)));
V  = zero_template;     % non-equilibrium factor
kd = zero_template;     % diss coef
for i = 1:M.num_elex_levels    % вот тут энергию e_i от нуля или от e_0?
 %  ZvibT=sum(exp(-M.ev_i{i}/(V_K*T)));
 %  ZvibU=sum(exp(M.ev_i{i}/(V_K*U(i))));
 ZvibT = sum( exp(-(M.ev_0(i) + M.ev_i{i}) / (V_K*T)) ); % vibr stat sum
 ZvibU = sum( exp( (M.ev_0(i) + M.ev_i{i}) / (V_K*U(i))) );
 divsum = divsum + M.s_e(i) * exp(M.e_E(i)/(V_K*U(i))) * ZvibU / ZvibT;
 V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
               exp((M.ev_0(i) + M.ev_i{i} + M.e_E(i))/V_K*(1/T + 1/U(i)));
end
V = Zel * V / divsum;
for i = 1:M.num_elex_levels
 if M.num_vibr_levels(i) > 1
  kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
      V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) ...
                                                               * kd_eq(i);
 else
  kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
                                                                kd_eq(i);
 end
end

end