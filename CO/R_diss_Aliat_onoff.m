function out = R_diss_Aliat_onoff(M, n_m, n_a1, n_a2, n_p, Coll, T, ...
                                            n0, ind_U, ind_Aliat)
% Универсальная функция расчёта диссоциации по Алиату для
% электронно-возбуждённых молекул.
% 22.10.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_U='3T';
   ind_Aliat='MT';
end

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
kd_eq = Coll.ArrA * T^Coll.ArrN*...
                            exp(-(M.diss_e + M.e_E)/(V_K*T)); % m^3/sec
% if ind_Arr==4
%     kd_eq = kd_eq*0+Coll.ArrA(ind_Arr);
% end

    % parameter of TM model
U=M.diss_e*0;
switch ind_U
    case 'D/6k'
        U = (M.diss_e+M.e_E)/(V_K*6);
    case '3T'
        U = U+3*T;
    case 'inf'
        U = U+Inf;
    case 'Savelev16'
        U = U+M.diss_e/V_K *(0.15+T/20000);
end
    % Aliat
% ind_Aliat stores data about Aliat model (1) or Marrone-Treanor model (0)
% to use.
if ind_Aliat=="Aliat"
 Zel=M.s_e(1:M.num_elex_levels)*exp(-M.e_E(1:M.num_elex_levels)'/(V_K*T));
else
 Zel=1;
end
divsum=0;
zero_template=zeros(1, sum(M.num_vibr_levels(1:M.num_elex_levels)));
V=zero_template;
Kdr=zero_template;
Theta_r = M.Be*V_H*V_C/V_K;
Z_rot = T./(M.sigma.*Theta_r);
for i=1:M.num_elex_levels       % вот тут энергию e_i от нуля или от e_0?
 ZvibT=sum(exp(-M.ev_i{i}/(V_K*T)));
 ZvibU=sum(exp(M.ev_i{i}/(V_K*U(i))));
 if ind_Aliat=="Aliat"
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
kd=zero_template;
for i=1:M.num_elex_levels
    kd(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
        V(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))...
                                                                *kd_eq(i);
end
kr= kd .* Kdr * n0;
RD = n_p * (n_a1*n_a2*kr-(n_m').*kd);
out = RD';
end