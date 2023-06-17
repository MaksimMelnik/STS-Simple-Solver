function out = R_exch_N2_O__NO_N(N2, NO, n_N2, n_O, n_NO, n_N, T)
% Функция расчёта обменной реакции N2(i)+O->NO(k)+N
% 28.05.2023

V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;

E=37484;
A = 3e-17^(T < 4000)*1.554e-23^(T >= 4000); % m^3/sec
b = 0^(T < 4000)*1.745^(T >= 4000);

%Park
% E=38370;
% N_a=6.02214076e23;
% A=6.4e17/N_a*1e-6;
% b=-1;

Theta_r_N2 = N2.Be(1)*V_H*V_C/V_K;     Z_rot_N2 = T./(N2.sigma.*Theta_r_N2);
Theta_r_NO = NO.Be(1)*V_H*V_C/V_K;     Z_rot_NO = T./(NO.sigma.*Theta_r_NO);

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку

%колебательные статсуммы
exp_n2 = exp(-(N2.ev_0(1) + N2.ev_i{1})/(V_K*T));
Zv_n2 = sum(exp_n2);

% приведенное Больцмановское распределение молекул
nn2i = exp_n2 / Zv_n2;


% разница энергий, K
dE_n2 = E - repmat((N2.ev_0(1) + N2.ev_i{1})'/V_K,1,NO.num_vibr_levels) + ...
    repmat((NO.ev_0(1)+NO.ev_i{1})/V_K,N2.num_vibr_levels,1);

EXP_n2 = exp(-dE_n2 .* heaviside(dE_n2) * (1/T));

kd_eq = A * T^b*exp(-E/T ); % m^3/sec

B_n2 = kd_eq / sum(EXP_n2 .* nn2i','all');

%скорости прямой реакции d-direct
kd = B_n2 * EXP_n2;

s_e_O_N=9/4; % s_e_O/s_e_N
m_O=2.6567628316576e-26;    m_N=2.32587E-26;

%отношение скорости обратной к скорости прямой
Kdr = (N2.mass*m_O/(NO.mass*m_N))^1.5*Z_rot_N2/Z_rot_NO*...
    exp((repmat(NO.ev_i{1}+NO.ev_0,N2.num_vibr_levels,1)-repmat((N2.ev_i{1}+N2.ev_0)',1,NO.num_vibr_levels))/(V_K*T))*...
    exp((N2.diss_e(1)-NO.diss_e(1))/V_K/T);
Kdr = N2.s_e(1)*s_e_O_N/NO.s_e(1) * Kdr;

%скорость обратной реакции r - reverse
kr= kd .* Kdr;

RExch1=n_NO'.*kr*n_N - n_N2.*kd*n_O;

%получаем матрицу 47 на 38, по строкам энергия N2, по столбцам NO
out = RExch1;
end