function out = R_diss_Aliat(M, n_m, n_a1, n_a2, n_p, Coll, T, n0, ...
    ind_Arr,ind_U)
% Универсальная функция расчёта диссоциации по Алиату для
% электронно-возбуждённых молекул.
% 22.10.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_Arr=1;
   ind_U=2;
end

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
kd_eq = Coll.ArrA(ind_Arr) * T^Coll.ArrN(ind_Arr)*...
    exp(-M.diss_e /(V_K*T)); % m^3/sec

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
    % Aliat
Zel=M.s_e*exp(-M.e_E'/(V_K*T));
divsum=0;
V=[];
Kdr=[];
Theta_r = M.Be*V_H*V_C/V_K;
Z_rot = T./(M.sigma.*Theta_r);
for i=1:M.num_elex_levels
    ZvibT=sum(exp(-M.e_i(i, 1:M.num_vibr_levels(i))/(V_K*T)));
    ZvibU=sum(exp(M.e_i(i, 1:M.num_vibr_levels(i))/(V_K*U)));
    divsum=divsum+M.s_e(i)*exp(M.e_E(i)/(V_K*U))*ZvibU/ZvibT;
    V=[V exp((M.e_i(i, 1:M.num_vibr_levels(i))+M.e_E(i))/V_K*(1/T+1/U))];
    Kdr=[Kdr Z_rot(i)*M.s_e(i)*...
        exp(-(M.e_i(i,1:M.num_vibr_levels(i))+M.e_E(i)-M.diss_e)/V_K/T)];
end
V=Zel*V/divsum;
Kdr=Kdr*(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2);
kd_new=V*kd_eq;

kd=kd_new;
Kdr2=0; %Kdr;
kr= kd .* Kdr2 * n0;
RD = n_p * (n_a1*n_a2*kr-(n_m').*kd);
out = RD';
end