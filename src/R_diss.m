function out = R_diss(M, n_m, n_a1, n_a2, n_p, Coll, T, n0, ind_Arr,ind_U)
% (T, ni_b, nm_b, nac_b, nao_b, naAr_b, modD,...
%     ind_Arr, n0)
% Универсальная функция расчёта диссоциации.
% Числовые плотности n_m и n_p получаем в безразмерном виде. 
% 10.10.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_Arr=1;
   ind_U=2;
end

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
kd_eq = Coll.ArrA(ind_Arr) * T^Coll.ArrN(ind_Arr)*...
    exp(-M.diss_e /(V_K*T));% m^3/sec

% if M.num_vibr_levels==1
%     Theta_r = M.Be*V_H*V_C/V_K;
%     Z_rot = T./(M.sigma.*Theta_r);
%     Z_vib=1;
%     Z_int=Z_rot*Z_vib;
%     Kdr2=(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2)...
%         *Z_int*exp(M.diss_e/V_K/T);
%     kr_eq= kd_eq .* Kdr2 * n0;
%     RD = n_p * (n_a1*n_a2*kr_eq-(n_m').*kd_eq);
% else
    ZvT = sum(exp(-M.e_i(1,:)'/(V_K*T)));
    % parameter of TM model
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
    ZvU = sum(exp(M.e_i(1,:)'/(V_K*U)));
    % non-equilibrium factor
    Z = ZvT / ZvU * exp((M.e_i(1,:)'+M.e_E(1))/V_K*(1/T + 1/U));
    % dis. rates
    kd = kd_eq' * Z'; % m^3/sec
    Theta_r = M.Be(1)*V_H*V_C/V_K;
    Z_rot = T./(M.sigma.*Theta_r);
    Kdr2=(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2)...
        *Z_rot *exp(-(M.e_i(1,:)+M.e_E(1)- M.diss_e)/V_K/T);
%     Kdr2=(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2)...
%         *Z_rot *exp(-(M.e_i- M.diss_e)/V_K/T);
    kr= kd .* Kdr2 * n0;
    RD = n_p * (n_a1*n_a2*kr-(n_m').*kd);
% end
out = RD';
end