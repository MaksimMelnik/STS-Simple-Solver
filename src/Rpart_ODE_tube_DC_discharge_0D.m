function out=Rpart_ODE_tube_DC_discharge_0D(~, y, kinetics)
% Right part function for ODE systems for the gas relaxation problem in DC
% discharge in a tube. The system is reduced to 0D [1].
% t is the time; y is the vector of macroparameters;
% kinetics is the big structure with all kinetics.
% 08.06.2023 Maksim Melnik
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.

    % constants
kb = 1.380649e-23;               % Boltzmann constant, J/K
N_a = 6.02214076e23;    % Avogadro constant
ne = y(end - 1);
%T_DN = y(end);          % dimentionless gas temperature T
T_DN = y(end - 1);
T = T_DN * kinetics.T0;                           % gas temperature T
Te = y(end);
Te = Te * kinetics.T0;

    % relaxation terms and energy flux
 y_rci = y(1:end-1);
[R, Q] = Rci(y_rci, kinetics);
R = R * kinetics.n0 * kinetics.t0;                % dimentionless 
Q = Q * kinetics.n0^2;                            % dimension value
Qe = 0;
if isKey(kinetics.reactions, 'free_e') %processes involving free electrons
 [Re, Qe] = Rci_e(y_rci, kinetics);
 Re = Re * kinetics.n0 * kinetics.t0;
 Qe = Qe * kinetics.n0^2;
 R = [R; 0] + Re;
end
Q_total = Q + Qe;
    % T equation
lambdaN2 = (1.717 + 0.084*T - 1.948e-5*T^2)/1e3;  % W / m / K
lambdaO2 = (1.056 + 0.087*T - 8.912e-6*T^2)/1e3;  % W / m / K
lambda = lambdaN2 * 0.8 + lambdaO2 * 0.2;         % W / m / K (kg*m/s3/K)
n_m = sum(y(1:end-1)) * kinetics.n0 / N_a;        % molar density, mol/m3
cp_N2 = c_p(kinetics.Ps{1}, T);
cp_O2 = c_p(kinetics.Ps{2}, T);
c_p_total = 0.8 * cp_N2 + 0.2 * cp_O2;      % molar heat capacity, J/mol/K
dT = (8*lambda*(kinetics.Tw - T)/kinetics.tube_R^2 + Q_total) ...
                        /(n_m*c_p_total) /kinetics.T0*kinetics.t0; % K/s
dTe = [];
if isKey(kinetics.reactions, 'free_e')
 dTe = (2 / 3 / kb ) * sum(Q) - R(end) * Te / ne;
end
% R_total = [R; 0] + Re;
% out = [R_total; dT];
out = [R; dT; dTe];
end