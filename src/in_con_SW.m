function [n1_DN, v1_DN, T1_DN] = in_con_SW(n0, v0, T0, rho0, f_m)
% The universal function for the calculation of macroparameters behind the
% shock wave (SW).
% n0 is the number density before a SW; v0 is the velocity before the SW;
% T is the gas temperature before the SW;
% rho0 is the density before a SW; f_m is the fraction of molecules.
% 17.06.2023

k = 1.380649e-23;
syms n1 v1 T1   % dim-less variables

C = v0^2 * rho0 / (n0 * k * T0);
S = [n1*v1 == 1, ...
    C*n1*v1^2 + n1*T1 == C + 1, ...
    (2.5+f_m)*T1 + 0.5*C*v1^2 == (2.5+f_m) + 0.5*C];
N = vpasolve(S, [n1,v1,T1], [5,100,0.2]);
Y2 = double([N.n1, N.v1, N.T1]);
ind = find(Y2(:, 3)>1);
n1_DN = Y2(ind(1), 1);
v1_DN = Y2(ind(1), 2);
T1_DN = Y2(ind(1), 3);
end