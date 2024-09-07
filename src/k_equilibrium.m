function kf = k_equilibrium(reaction, T)
% Calculation of equilibrium rate coefficient k (m3/s).
% reaction is the reaction variable with required parameters; 
% T is the gas temperature. 
% 06.09.2024 by Maksim Melnik
k = 1.380649e-23;         % Boltzmann constant, J/K
 switch reaction.type    % choosing reaction type
  case "const"
    kf = reaction.A;
  case "ATn"
    kf = reaction.A * T ^ reaction.n;
  case "Arrhenius"
    kf = reaction.A * T ^ reaction.n * exp(- reaction.E / k / T);
  case "A(T/d_T)^n"
    kf = reaction.A * (T / reaction.d_T) ^ reaction.n * ...
                                                exp(- reaction.E / k / T);
  case "Starik_test_on_T"
    kf = reaction.A(T) * T ^ reaction.n(T) * exp(- reaction.E(T) / k /T);
  otherwise
    error("Equilibrium rate coefficients of this type are still not " + ...
            "implemented: " + reaction.type)
 end
end