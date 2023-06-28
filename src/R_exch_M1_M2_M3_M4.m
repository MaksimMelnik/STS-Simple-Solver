function out = R_exch_M1_M2_M3_M4(M1, M2, M3, M4, n_M1, n_M2, n_M3, n_M4, T, coll, VDOP)
% Функция расчёта обменной реакции M1+M2->M3+M4
%M1 и M3 молекулы с колебательными уровнями 
%VDOP - Vibrational deactivation of product. 1 - учитываем активацию
%продукта реакции, 0 - не учитываем.
% 28.06.2023

V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34; %константы

%Параметры Аррениуса в переменной Coll: ArrA, ArrN, ArrE
%(kd_eq=A*T^N*exp(-E/T))

Theta_r_M1 = M1.Be(1)*V_H*V_C/V_K;     Z_rot_M1 = T./(M1.sigma.*Theta_r_M1);
Theta_r_M3 = M3.Be(1)*V_H*V_C/V_K;     Z_rot_M3 = T./(M3.sigma.*Theta_r_M3);

%колебательные статсуммы
exp_M1 = exp(-(M1.ev_0(1) + M1.ev_i{1}(1:M1.num_vibr_levels(1)))/(V_K*T));
Zv_M1 = sum(exp_M1);

% приведенное Больцмановское распределение молекул
nM1i = exp_M1 / Zv_M1;


% разница энергий, K
dE_M1 = coll.ArrE - repmat((M1.ev_0(1) + M1.ev_i{1}(1:M1.num_vibr_levels(1)))'/V_K,1,M3.num_vibr_levels(1)) + ...
    repmat((M3.ev_0(1)+VDOP*M3.ev_i{1}(1:M3.num_vibr_levels(1)))/V_K, M1.num_vibr_levels(1), 1);

EXP_M1 = exp(-dE_M1 .* heaviside(dE_M1) * (1/T));

kd_eq = coll.ArrA * T^coll.ArrN*exp(-coll.ArrE/T ); % m^3/sec

B_M1 = kd_eq / sum(EXP_M1 .* nM1i','all');

%скорости прямой реакции d-direct
kd = B_M1 * EXP_M1;

%отношение скорости обратной к скорости прямой
Kdr = (M1.mass*M2.mass/(M3.mass*M4.mass))^1.5*Z_rot_M1/Z_rot_M3*...
exp((repmat(M3.ev_i{1}(1:M3.num_vibr_levels(1))*VDOP+M3.ev_0,M1.num_vibr_levels(1),1) ...
-repmat((M1.ev_i{1}(1:M1.num_vibr_levels(1))+M1.ev_0)', 1 , M3.num_vibr_levels(1) ))/(V_K*T))*...
exp((M1.diss_e(1)-M3.diss_e(1))/V_K/T);
Kdr = M1.s_e(1)*M2.s_e(1)/(M3.s_e(1)*M4.s_e(1)) * Kdr;

%скорость обратной реакции r - reverse
kr= kd .* Kdr;
RExch1=n_M3'.*kr*n_M4 - n_M1.*kd*n_M2;

%получаем матрицу M1.num_vibr_levels на M3.num_vibr_levels, по строкам
%энергия M1, по столбцам M3
out = RExch1;
end