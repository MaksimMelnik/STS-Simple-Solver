function out=k_diss_old(M, Coll, T, n0, ind_Arr, ind_U)
% коэффициент диссоциации по Алиату
% 24.10.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<5
   ind_Arr=1;
   ind_U=2;
end
kd_eq = Coll.ArrA(ind_Arr) * T.^Coll.ArrN(ind_Arr).*...
    exp(-M.diss_e./(V_K*T)); % m^3/sec
switch ind_U
    case 2
        U = M.diss_e/(V_K*6);
    case 3
        U = 3*T;
    case 4
        U = Inf;
    case 5
        U = M.diss_e/V_K *(0.15+T/20000);
end
ZvT = sum(exp(-M.e_i(1,:)'./(V_K*T)),1);
    % Aliat
ZvU = sum(exp(M.e_i(1,:)'./(V_K*U)),1);
    % non-equilibrium factor
Z = ZvT ./ ZvU .* exp((M.e_0(1)+M.e_i(1,:)')/V_K.*(1./T + 1./U));
    % dis. rates
kd = kd_eq .* Z; % m^3/sec
Theta_r = M.Be(1)*V_H*V_C/V_K;
Z_rot = T./(M.sigma.*Theta_r);
Kdr2=(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3.*(2*pi*V_K*T).^(-3/2)...
    .*Z_rot .*exp(-(M.e_i(1,:)'+M.e_E(1)- M.diss_e)/V_K./T);
kr= kd .* Kdr2;
out=Z;
end