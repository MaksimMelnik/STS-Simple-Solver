function out=Rpart_0D(t, y, kinetics)%#ok<INUSL>


% constants
k = 1.380649e-23;       % Boltzmann constant, J/K
T_DN = y(end);          % dimentionless gas temperature T

% relaxation terms
[R, ~] = Rci(y, kinetics);
R = R * kinetics.n0 * kinetics.t0;
    % number densities equations/
M = R;
% energy equation
nm = 0; % number density of molecules
na = 0; 
Rci_term = 0;
E_term = 0;
for ind = 1:kinetics.num_Ps
    if kinetics.Ps{ind}.fr_deg_c > 3
        e_i = [];
        for ind_e = 1:kinetics.Ps{ind}.num_elex_levels
            e_i = [e_i, kinetics.Ps{ind}.e_E(ind_e) + kinetics.Ps{ind}.ev_i{ind_e} + kinetics.Ps{ind}.ev_0(ind_e)];
        end
        E_term = E_term + sum(R(kinetics.index{ind}) .* (e_i + kinetics.Ps{ind}.form_e)) / k / kinetics.T0;
        Rci_term = Rci_term + 2.5 * sum(R(kinetics.Ps{ind}));
        nm = nm + sum(y(kinetics.index{ind}));
    else
        E_term = E_term + R(kinetics.Ps{ind}) * (kinetics.Ps{ind}.form_e + kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels)) / k / kinetics.T0;
        Rci_term = Rci_term + 1.5 * sum(R(kinetics.Ps{ind}));
        na = na + sum(y(kinetics.index{ind}));
    end
end
M(kinetics.num_eq + 1) = -2 * (T_DN * Rci_term + E_term) / (5 * nm + 3 * na);
out = M;
end
