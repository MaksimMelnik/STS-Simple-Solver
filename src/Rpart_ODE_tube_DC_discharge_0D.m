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
% ne = y(end - 1);
ne = y(end - 2);         % dimentionless
%T_DN = y(end);          % dimentionless gas temperature T
T_DN = y(end - 1);
T = T_DN * kinetics.T0;                           % gas temperature T
Te = y(end)* kinetics.T0;
% Te = Te * kinetics.T0;

    % relaxation terms and energy flux
 y_rci = y(1:end-1);
[Rh, Q] = Rci(y_rci, kinetics);                 % heavy particles
R = Rh * kinetics.n0 * kinetics.t0;             % dimentionless 
Q = Q * kinetics.n0^2;                          % dimension value, kg/m/s3
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
                /(n_m*c_p_total) /kinetics.T0*kinetics.t0; % dimentionless
dTe = [];
if isKey(kinetics.reactions, 'free_e')
 me = 9.1094e-31;
 fr = 0;
 for i = [1, 2]
    mi = kinetics.Ps{i}.mass;
    mred = mi*me/(mi+me);
    ni = sum(y(kinetics.index{i}));
    r = kinetics.Ps{i}.diameter / 2;
    z = sqrt(8*pi*kb*T/mred)*r^2;    % m3/s
    fr = fr + z*ni/mi*kinetics.n0;      % 1/kg/s
 end
 dTe = (2 / 3 / kb) * Qe / ne / kinetics.T0 / kinetics.n0 * kinetics.t0...
             - R(end) * Te / ne / kinetics.T0 ...
             - 3*(Te-T)*fr*ne*me*kinetics.t0/kinetics.T0;  % dimentionless
 % dTe = 0;
end
% R_total = [R; 0] + Re;
% out = [R_total; dT];
out = [R; dT; dTe];
end