function check_CL_0D(L0, Y, kinetics, ind_ref)
% The function aimed to check if conservation laws are respected in the 
% 0D relaxation problem.
% L0 is the initial values; Y is the evolution of macroparameters; 
% kinetics is the kinetics container; ind_ref is the indicator for
% the reference value (0 --- L0 values, 1 --- the first values in Y).
% 02.06.2023 Maksim Melnik

    % constants
k = 1.380649e-23;             % Boltzmann constant, J/K
rho = 0;
E = 0;
for ind=1:kinetics.num_Ps
    % rho
 rho=rho + sum(Y(:, kinetics.index{ind}), 2) * kinetics.Ps{ind}.mass;
    % E
 if kinetics.Ps{ind}.fr_deg_c > 3
  e_i = [];
  for ind_e = 1:kinetics.Ps{ind}.num_elex_levels
   e_i = [e_i, kinetics.Ps{ind}.e_E(ind_e) + ...
            kinetics.Ps{ind}.ev_i{ind_e} + kinetics.Ps{ind}.ev_0(ind_e)];
  end
  E = E + sum(Y(:, kinetics.index{ind}) .* ...
                    (e_i + kinetics.Ps{ind}.form_e + 2.5*k*Y(:, end)), 2);
 else
  E = E + sum(Y(:, kinetics.index{ind}).* ...
       (kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels)+...
                            kinetics.Ps{ind}.form_e+1.5*k*Y(:,end)), 2);
 end
end
E0 = L0(2);
rho0 = L0(1);
if ind_ref == 1
    rho0=rho(1);
    E0 = E(1);
end
disp(['rho max error ' num2str(max(abs( (rho-rho0)/rho0 )))])
disp(['E cons max error ' num2str(max(abs( (E-E0)/E0 )))])

maxerror = 9e-5;
if max(abs( (rho-rho0)/rho0 )) > maxerror 
    warning("rho max error")
end
if max(abs( (E-E0)/E0 )) > maxerror 
    warning("E cons max error")
end
end