function out = R_exch_CO_C__C2_O(C2, CO, n_m, n_C, n_O, n_C2, Coll, T, ...
                                            ...n0, 
                                            ind_Arr)%, ind_U, ind_Aliat)
% Функция расчёта обменной реакции CO+C->C2+O
% 11.07.2022
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_Arr=1;
%    ind_U=2;
%    ind_Aliat=0;
end


% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
Theta_r_CO = CO.Be*V_H*V_C/V_K;     Z_rot_CO = T./(CO.sigma.*Theta_r_CO);
Theta_r_C2 = C2.Be*V_H*V_C/V_K;     Z_rot_C2 = T./(C2.sigma.*Theta_r_C2);
Kdr=zeros(1, sum(CO.num_vibr_levels));
s_e_C_O=1/9; % s_e_C/s_e_O
m_O=2.6567628316576e-26;    m_C=1.994473440944e-26;
if (ind_Arr==4)||(ind_Arr==7)
%     C2+O->CO+C
    kd_rec = Coll.ArrA(ind_Arr);
%     не, тут вообще нужно реакцию наоборот вывернуть
   for i=1:CO.num_elex_levels
  Kdr(1+sum(CO.num_vibr_levels(1:i-1)):sum(CO.num_vibr_levels(1:i)))=...
        exp(-(CO.diss_e(i)-CO.e_E(i)-CO.ev_0(i)-CO.ev_i{i})/V_K/T)...
        /Z_rot_CO(i)/CO.s_e(i);
   end
    Kdr=Kdr/s_e_C_O*C2.s_e*(C2.mass*m_O/CO.mass/m_C)^1.5*Z_rot_C2*...
        exp((C2.diss_e)/V_K/T);
    kd=kd_rec.*Kdr;
    kr=kd_rec;
else
 for i=1:CO.num_elex_levels
    Kdr(1+sum(CO.num_vibr_levels(1:i-1)):sum(CO.num_vibr_levels(1:i)))=...
        Z_rot_CO(i)*CO.s_e(i)... 
        *exp((CO.diss_e(i)-CO.e_E(i)-CO.ev_0(i)-CO.ev_i{i})/V_K/T);
 end
    Kdr=Kdr*s_e_C_O/C2.s_e*(CO.mass*m_C/C2.mass/m_O)^1.5/Z_rot_C2*...
                                                exp((-C2.diss_e)/V_K/T);
    kd_eq = Coll.ArrA(ind_Arr) * T^Coll.ArrN(ind_Arr)*...
                                exp(-(C2.CO_C_C2_O_e)/(V_K*T)); % m^3/sec
    kr= kd_eq .* Kdr;
    kd=kd_eq;
end
% disp(Kdr)
RExch=n_C2*n_O*kr -n_m'*n_C.*kd;
out = RExch';
%     kd_eq = Coll.ArrA(1) * T^Coll.ArrN(1)*...
%                                 exp(-(C2.CO_C_C2_O_e)/(V_K*T)); % m^3/sec
%     kr= kd_eq .* Kdr;
%     kd=kd_eq;
%     RExch2=n_C2*n_O*kr -n_m'*n_C.*kd;
%     disp([num2str(RExch(1)) ' ' num2str(RExch2(1))])
end