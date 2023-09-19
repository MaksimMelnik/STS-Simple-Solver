function kd = k_diss_Marrone_Treanor(M, T, kd_eq, U)
% Dissociation coefficient using the Marrone-Treanor model.
% kd is the vector of rate coefficients.
% M is the dissociating molecule; T is the gas temperature;
% kd_eq is the equlibrium rate coefficient;
% U is the non-equilibrium parameter.
% 19.09.2023 Maksim Melnik

k = 1.380649e-23;   % Boltzmann constant
zero_template = zeros(1, sum(M.num_vibr_levels(1:M.num_elex_levels)));
V  = zero_template; % non-equilibrium factor
kd = zero_template; % diss coef
for i = 1:M.num_elex_levels    % вот тут энергию e_i от нуля или от e_0?
 if M.num_vibr_levels(i) > 1
  %  ZvibT=sum(exp(-M.ev_i{i}/(k*T)));
  %  ZvibU=sum(exp(M.ev_i{i}/(k*U(i))));
  ZvibT = sum( exp(-(M.ev_0(i) + M.ev_i{i})/(k*T)) ); %vibrational statsum
  ZvibU = sum( exp( (M.ev_0(i) + M.ev_i{i})/(k*U(i))) );
  % V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=... old
  %                exp((M.ev_i{i}+M.e_E(i))/k*(1/T+1/U(i)))*ZvibT/ZvibU;
  V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...new
        exp((M.ev_0(i) + M.ev_i{i})/k*(1/T + 1/U(i))) * ZvibT / ZvibU;
  kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
      V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) ...
                                                               * kd_eq(i);
 else
  kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i))) = ...
                                                                kd_eq(i);
 end
end

end