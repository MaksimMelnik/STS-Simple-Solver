function out = R_diss_Aliat_onoff(M, n_m, n_a1, n_a2, n_p, Coll, T, ...
                                            n0, ind_Arr, ind_U, ind_Aliat)
% Универсальная функция расчёта диссоциации по Алиату для
% электронно-возбуждённых молекул.
% 22.10.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_Arr=1;
   ind_U=2;
   ind_Aliat=0;
end

% disp(Coll.ArrA)
% pause

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
kd_eq = Coll.ArrA(ind_Arr) * T^Coll.ArrN(ind_Arr)*...
                            exp(-(M.diss_e + M.e_E)/(V_K*T)); % m^3/sec
if ind_Arr==4
    kd_eq = kd_eq*0+Coll.ArrA(ind_Arr);
end

    % parameter of TM model
U=M.diss_e*0;
switch ind_U
    case 2
        U = (M.diss_e+M.e_E)/(V_K*6);
    case 3
        U = U+3*T;
    case 4
        U = U+Inf;
    case 5
        U = U+M.diss_e/V_K *(0.15+T/20000);
end
    % Aliat
% ind_Aliat stores data about Aliat model (1) or Marrone-Treanor model (0)
% to use.
if ind_Aliat
 Zel=M.s_e(1:M.num_elex_levels)*exp(-M.e_E(1:M.num_elex_levels)'/(V_K*T));
else
 Zel=1;
end
divsum=0;
V=zeros(1, sum(M.num_vibr_levels));
Kdr=zeros(1, sum(M.num_vibr_levels));
Theta_r = M.Be*V_H*V_C/V_K;
Z_rot = T./(M.sigma.*Theta_r);
for i=1:M.num_elex_levels       % вот тут энергию e_i от нуля или от e_0?
 ZvibT=sum(exp(-M.ev_i{i}/(V_K*T)));
 ZvibU=sum(exp(M.ev_i{i}/(V_K*U(i))));
 if ind_Aliat
  divsum=divsum+M.s_e(i)*exp(M.e_E(i)/(V_K*U(i)))*ZvibU/ZvibT;
  V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
                exp((M.ev_i{i}+M.e_E(i))/V_K*(1/T+1/U(i)));
 else
  V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
                exp((M.ev_i{i}+M.e_E(i))/V_K*(1/T+1/U(i)))*ZvibT/ZvibU;
  divsum=1;
 end
 Kdr(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
     Z_rot(i)*M.s_e(i)/M.mltpl_atoms_s_e... 
     *exp(-(M.e_E(i)+M.ev_0(i)+M.ev_i{i}...
                                   +M.form_e-M.form_e_atoms_sum)/V_K/T);
end
V=Zel*V/divsum;
Kdr=Kdr*(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2);
kd=zeros(1, sum(M.num_vibr_levels));
for i=1:M.num_elex_levels
    kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
        V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))...
                                                                *kd_eq(i);
end
kr= kd .* Kdr * n0;
% kr= 0;
RD = n_p * (n_a1*n_a2*kr-(n_m').*kd);
out = RD';
%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('kk')
% disp(kr)
% disp(kd_eq)
% disp(abs(diff)>0)
% disp(kd_eq)
% disp('MT')
% disp(kd_eq(1))
% warning('R_diss_Aliat_onoff returns kr')
% out = kr/n0;
end