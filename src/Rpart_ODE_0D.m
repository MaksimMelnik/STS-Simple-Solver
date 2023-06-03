function out=Rpart_ODE_0D(t, y, kinetics)%#ok<INUSL>
% Right part function for ODE systems for 0D relaxation problem.
% t is the time; y is the vector of macroparameters;
% kinetics is the big structure with all kinetics.
% 31.05.2023 Maksim Melnik

    % constants
k=1.380649e-23;     % Boltzmann constant, J/K
T_DN=y(end);        % dimentionless gas temperature T

    % relaxation terms
R=Rci(y, kinetics);
    % number densities equations
M = diag([ones(1, kinetics.num_eq) 0]);
    % energy equation
nm=0;   % number density of molecules
for ind = 1:kinetics.num_Ps
 if kinetics.Ps{ind}.fr_deg_c > 3
  e_i = [];
  for ind_e = 1:kinetics.Ps{ind}.num_elex_levels
   e_i = [e_i, kinetics.Ps{ind}.e_E(ind_e) + ...
            kinetics.Ps{ind}.ev_i{ind_e} + kinetics.Ps{ind}.ev_0(ind_e)];
  end
  M(end, kinetics.index{ind}) = ...
                (e_i + kinetics.Ps{ind}.form_e)/k/kinetics.T0 + 2.5*T_DN;
  nm = nm + sum(y(kinetics.index{ind}));
 else
  M(end, kinetics.index{ind}) = (kinetics.Ps{ind}.form_e ...
      + kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels))...
                                                /k/kinetics.T0 + 1.5*T_DN;
 end
end
M(end, end) = 1.5*sum(y(1:end-1)) + nm;
Msp = sparse(M);
R = [R; 0] * kinetics.n0 * kinetics.t0;
out = Msp^(-1) * R;

% out=y*0;

end