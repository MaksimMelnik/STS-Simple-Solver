function out = R_exch_N2_O__NO_N(N2, NO, n_N2, n_O, n_NO, n_N, Coll, T)
% Функция расчёта обменной реакции N2+O->NO+N
% 28.05.2023
%написана по полной аналогии реакции NO + O -> O2 + N
%заменой NO на N2 и O2 на NO

V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

E=38370;

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку

Theta_r_N2 = N2.Be(1)*V_H*V_C/V_K;     Z_rot_N2 = T./(N2.sigma.*Theta_r_N2);
Theta_r_NO = NO.Be(1)*V_H*V_C/V_K;     Z_rot_NO = T./(NO.sigma.*Theta_r_NO);
Kdr=zeros(sum(NO.num_vibr_levels), sum(N2.num_vibr_levels));
s_e_O_N=9/4; % s_e_O/s_e_N
m_O=2.6567628316576e-26;    m_N=2.32587E-26;

i=1;
for j=1:NO.num_vibr_levels
Kdr(j, 1+sum(N2.num_vibr_levels(1:i-1)):sum(N2.num_vibr_levels(1:i)))=...
        Z_rot_N2*N2.s_e(1)... 
        *exp((N2.diss_e(1) - NO.diss_e(1) + NO.ev_0(1)  + NO.ev_i{1}(j) -...
        N2.ev_0(1) - N2.ev_i{1})/V_K/T)...
        *s_e_O_N/NO.s_e(1)*(N2.mass*m_O/NO.mass/m_N)^1.5/Z_rot_NO;
end
    kd_eq = Coll.ArrA(1) * T^Coll.ArrN(1)*...
                                exp(-E/T ); % m^3/sec
    kr= kd_eq .* Kdr;
  
    kd=kd_eq;

RExch1=n_NO.*kr*n_N - n_N2'.*n_O*kd;
%получаем матрицу 47 на 47, по строкам энергия NO, по столбцам N2
out = RExch1;
end