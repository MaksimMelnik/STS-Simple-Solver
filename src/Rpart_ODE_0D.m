function out = Rpart_ODE_0D(~, y, kinetics)
% Right part function for ODE systems for 0D relaxation problem.
% y is the vector of macroparameters;
% kinetics is the big structure with all kinetics.
% 21.07.2023 Shaikhutdinova Asya

    % constants
k = 1.380649e-23;       % Boltzmann constant, J/K
T_DN = y(end);          % dimentionless gas temperature T

    % relaxation terms
[R, ~] = Rci(y, kinetics);
R = R * kinetics.n0 * kinetics.t0;

    % number densities equations
out = [R; 0];

    % energy equation
nm = 0;                 % number density of molecules
Rci_term = 0;
E_term = 0;
for ind = 1:kinetics.num_Ps
 if kinetics.Ps{ind}.fr_deg_c > 3
  e_i = [];
  for ind_e = 1:kinetics.Ps{ind}.num_elex_levels
   e_i = [e_i, kinetics.Ps{ind}.e_E(ind_e) + ...
            kinetics.Ps{ind}.ev_i{ind_e} + kinetics.Ps{ind}.ev_0(ind_e)];
  end
  E_term = E_term + ...
        sum((R(kinetics.index{ind}))' .* (e_i + kinetics.Ps{ind}.form_e));
  Rci_term = Rci_term + sum(R(kinetics.index{ind}));
  nm = nm + sum(y(kinetics.index{ind}));
 else
  E_term = E_term + ...
      sum(R(kinetics.index{ind}) .* (kinetics.Ps{ind}.form_e + ...
            kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels)));
 end
end
Rci_term = Rci_term + 1.5 * sum(R);
n_all = sum(y(1:end-1));
out(end) = - (T_DN * Rci_term + E_term / k / kinetics.T0) / ...
                                                    (nm + 1.5 * n_all);
end