function [R, Q] = R_ions_wall(M1, M2, n_a, T, kinetics)
        % e + O2+(X) + wall -> O2(X)
% to finish

switch M1.name
    case "O2+"
        mu_i = 2.24 / 1e4;  % m^2/V/s, for O2+
    case "N2+"
        mu_i = 1.87 / 1e4;
    otherwise
            
end
mu_i = mu_i * 1e3;                      % temporal correction
eV_to_J = 1.60218e-19;
e_charge = 1.602176634e-19;             % in C
Lambda_i = kinetics.tube_R/2.405;       % m
u_k = 2.59294906071072e+00 * eV_to_J;   % J
D_ai = mu_i * u_k / e_charge;
R = D_ai / Lambda_i^2 * n_a;
Q = (M1.form_e - M2.form_e) * R;
end