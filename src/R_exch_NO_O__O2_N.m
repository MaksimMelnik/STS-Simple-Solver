function out = R_exch_NO_O__O2_N(NO, O2, n_NO, n_O, n_O2, n_N,Coll, T)
% Функция расчёта обменной реакции NO+O->O2+N
% 28.05.2023

V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

E=19450;

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку

Theta_r_NO = NO.Be(1)*V_H*V_C/V_K;     Z_rot_NO = T./(NO.sigma.*Theta_r_NO);
Theta_r_O2 = O2.Be(1)*V_H*V_C/V_K;     Z_rot_O2 = T./(O2.sigma.*Theta_r_O2);
Kdr=zeros(sum(O2.num_vibr_levels), sum(NO.num_vibr_levels));
s_e_O_N=9/4; % s_e_O/s_e_N
m_O=2.6567628316576e-26;    m_N=2.32587E-26;

i=1;
for j=1:O2.num_vibr_levels
Kdr(j, 1+sum(NO.num_vibr_levels(1:i-1)):sum(NO.num_vibr_levels(1:i)))=...
        Z_rot_NO*NO.s_e(1)... 
        *exp((NO.diss_e(1) - O2.diss_e(1) + O2.ev_0(1)  + O2.ev_i{1}(j) -...
        NO.ev_0(1) - NO.ev_i{1})/V_K/T)...
        *s_e_O_N/O2.s_e(1)*(NO.mass*m_O/O2.mass/m_N)^1.5/Z_rot_O2;
end
    kd_eq = Coll.ArrA(1) * T^Coll.ArrN(1)*...
                                exp(-E/T ); % m^3/sec
    kr= kd_eq .* Kdr;
  
    kd=kd_eq;
% disp(n_O2);
% disp(kr);
% disp(n_O2.*kr*n_N);
 %disp(n_NO'.*n_O*kd);
RExch1=n_O2.*kr*n_N - n_NO'.*n_O*kd;
%disp(RExch1);
out = RExch1;
end