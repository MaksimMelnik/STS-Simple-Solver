function [RExch1, Q] = R_exch_O2_N__NO_O(O2, NO, n_O2, n_N, n_NO, n_O, T)
% Функция расчёта обменной реакции O2+N->NO+O
% 28.05.2023

V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34; %константы

%Значения взяты из статьи
A = 4e-16^(T < 4000)*3.206e-23^(T >= 4000); % m^3/sec

b = (-0.39)^(T < 4000)*1.58^(T >= 4000);

E =  1449; % K

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку

Theta_r_NO = NO.Be(1)*V_H*V_C/V_K;     Z_rot_NO = T./(NO.sigma.*Theta_r_NO);
Theta_r_O2 = O2.Be(1)*V_H*V_C/V_K;     Z_rot_O2 = T./(O2.sigma.*Theta_r_O2);

%колебательные статсуммы
exp_o2 = exp(-(O2.ev_0(1) + O2.ev_i{1})/(V_K*T));
Zv_o2 = sum(exp_o2);

% приведенное Больцмановское распределение молекул
no2i = exp_o2 / Zv_o2;


% разница энергий, K
dE_o2 = E - repmat((O2.ev_0(1) + O2.ev_i{1})'/V_K,1,NO.num_vibr_levels) + ...
    repmat((NO.ev_0(1)+NO.ev_i{1})/V_K,O2.num_vibr_levels,1);

EXP_o2 = exp(-dE_o2 .* heaviside(dE_o2) * (1/T));

kd_eq = A * T^b*exp(-E/T ); % m^3/sec

B_o2 = kd_eq / sum(EXP_o2 .* no2i','all');

%скорости прямой реакции d-direct
kd = B_o2 * EXP_o2;

s_e_N_O=4/9; % s_e_N/s_e_O
m_O=2.6567628316576e-26;    m_N=2.32587E-26;

%отношение скорости обратной к скорости прямой
dE = (repmat(NO.ev_i{1} + NO.ev_0, O2.num_vibr_levels, 1) - ...
            repmat((O2.ev_i{1} + O2.ev_0)', 1, NO.num_vibr_levels)) + ...
    (O2.diss_e(1)-NO.diss_e(1));
Kdr = (O2.mass*m_N/(NO.mass*m_O))^1.5 * Z_rot_O2 / Z_rot_NO * ...
    exp( dE / (V_K*T));
Kdr = O2.s_e(1)*s_e_N_O/NO.s_e(1) * Kdr;

%скорость обратной реакции r - reverse
kr= kd .* Kdr;
RExch1=n_NO'.*kr*n_O - n_O2.*kd*n_N;

Q = sum(RExch1 .* dE, 'all');

%получаем матрицу 36 на 38, по строкам энергия O2, по столбцам NO
% out = RExch1;
end