function out=Rpart_ODE_tube_DC_discharge_0D(t, y, kinetics)%#ok<INUSL>
% Right part function for ODE systems for the gas relaxation problem in DC
% discharge in a tube. The system is reduced to 0D [1].
% t is the time; y is the vector of macroparameters;
% kinetics is the big structure with all kinetics.
% 08.06.2023 Maksim Melnik
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.

    % constants
N_a = 6.02214076e23;    % Avogadro constant
T_DN = y(end);          % dimentionless gas temperature T
T = T_DN * kinetics.T0;                           % gas temperature T

    % relaxation terms and energy flux
[R, Q] = Rci(y, kinetics);
R = R * kinetics.n0 * kinetics.t0;                % dimentionless 
Q = Q * kinetics.n0^2;                            % dimension value
    % T equation
lambdaN2 = (1.717 + 0.084*T - 1.948e-5*T^2)/1e3;  % W / m / K
lambdaO2 = (1.056 + 0.087*T - 8.912e-6*T^2)/1e3;  % W / m / K
lambda = lambdaN2 * 0.8 + lambdaO2 * 0.2;         % W / m / K (kg*m/s3/K)
n_m = sum(y(1:end-1)) * kinetics.n0 / N_a;        % molar density, mol/m3
cp_N2 = 29.1 + 2494.2/553.4/sqrt(pi/2) * exp(-2 * ((T-1047.4)/553.4)^2);
cp_O2 = 28.2 + 6456.2/788.3/sqrt(pi/2) * exp(-2 * ((T-1006.9)/788.3)^2);
c_p = 0.8 * cp_N2 + 0.2 * cp_O2;            % molar heat capacity, J/mol/K
dT = (8*lambda*(kinetics.Tw - T)/kinetics.tube_R^2 + Q) / (n_m*c_p) ...
                                        /kinetics.T0*kinetics.t0; % K/s
out = [R; dT];
end