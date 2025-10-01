% Old ptau (pressure * vibrational relaxation time) plot for O2-O
% original 22.05.2020
% 15.09.2025 by Maksim Melnik
%% Content:
% 1) p*tau for O2-O
% 2) p*tau for NO-NO
% 3) p*tau for CO-CO

addpath('../src/')
load("../data/particles.mat", "O2", "O", "NO", "CO")
h = 6.6261*10^(-34);
k = 1.3807e-23;
c = 2.99*10^10;

%% 1) p*tau for O2-O
% данные из каппы
gg=[
1	0	300		3.20193e-23	2.87047e-18	3.5954e-18	3.91896e-18	9.94356e-18	
1	0	400		1.33909e-22	3.43365e-18	3.97556e-18	4.0301e-18	1.49798e-17	
1	0	500		4.05261e-22	3.97192e-18	4.3234e-18	4.13799e-18	2.00568e-17	
1	0	600		1.01689e-21	4.49768e-18	4.6521e-18	4.23985e-18	2.50641e-17	
1	0	700		2.19764e-21	5.01825e-18	4.96899e-18	4.33553e-18	2.99517e-17	
1	0	800		4.24741e-21	5.53839e-18	5.2786e-18	4.42549e-18	3.46979e-17	
1	0	900		7.53162e-21	6.06139e-18	5.58394e-18	4.51033e-18	3.9295e-17	
1	0	1000	1.24761e-20	6.58964e-18	5.8871e-18	4.59063e-18	4.37429e-17	
1	0	1500	8.04131e-20	9.36179e-18	7.41243e-18	4.93934e-18	6.39172e-17	
1	0	2000	2.77125e-19	1.24292e-17	9.00977e-18	5.22655e-18	8.12254e-17	
1	0	2500	6.87256e-19	1.58487e-17	1.07144e-17	5.4734e-18	9.63427e-17	
1	0	3000	1.39548e-18	1.96496e-17	1.25405e-17	5.69149e-18	1.09755e-16	
1	0	3500	2.48058e-18	2.38473e-17	1.44934e-17	5.88791e-18	1.21806e-16	
1	0	4000	4.01289e-18	2.84491e-17	1.65741e-17	6.06731e-18	1.32747e-16	
1	0	4500	6.05371e-18	3.34559e-17	1.87805e-17	6.23294e-18	1.42766e-16	
1	0	5000	8.65554e-18	3.88644e-17	2.1109e-17	6.38715e-18	1.52006e-16	
1	0	5500	1.18627e-17	4.46678e-17	2.35549e-17	6.53172e-18	1.60578e-16	
1	0	6000	1.57123e-17	5.08564e-17	2.61125e-17	6.66802e-18	1.68573e-16	
1	0	6500	2.02349e-17	5.74182e-17	2.87756e-17	6.79712e-18	1.76061e-16	
1	0	7000	2.54555e-17	6.43395e-17	3.15377e-17	6.91991e-18	1.83104e-16	
1	0	7500	3.1394e-17	7.16051e-17	3.4392e-17	7.0371e-18	1.89751e-16	
1	0	8000	3.80663e-17	7.9199e-17	3.73316e-17	7.14929e-18	1.96042e-16	
1	0	8500	4.54845e-17	8.71041e-17	4.03496e-17	7.25696e-18	2.02014e-16	
1	0	9000	5.36575e-17	9.5303e-17	4.34392e-17	7.36056e-18	2.07697e-16	
1	0	9500	6.25914e-17	1.03778e-16	4.65937e-17	7.46044e-18	2.13117e-16	
1	0	10000	7.22901e-17	1.1251e-16	4.98063e-17	7.55691e-18	2.18296e-16	
1	0	11000	9.39877e-17	1.30678e-16	5.63807e-17	7.74071e-18	2.28012e-16	
1	0	12000	1.18745e-16	1.49662e-16	6.31135e-17	7.91378e-18	2.36976e-16	
1	0	13000	1.46538e-16	1.69326e-16	6.99597e-17	8.07755e-18	2.45294e-16	
1	0	14000	1.77325e-16	1.8954e-16	7.68779e-17	8.23317e-18	2.53051e-16	
1	0	15000	2.11056e-16	2.1018e-16	8.38312e-17	8.38159e-18	2.60314e-16	
1	0	16000	2.47671e-16	2.31133e-16	9.07863e-17	8.52357e-18	2.67142e-16	
1	0	17000	2.87106e-16	2.52294e-16	9.77138e-17	8.65977e-18	2.73581e-16	
1	0	18000	3.29295e-16	2.7357e-16	1.04588e-16	8.79073e-18	2.79672e-16	
1	0	19000	3.74167e-16	2.94873e-16	1.11387e-16	8.91694e-18	2.85449e-16	
1	0	20000	4.21652e-16	3.16127e-16	1.1809e-16	9.0388e-18	2.90941e-16	
];

ptau_Torres = [
% # --O2-O:Minnesota PES# -------combined-----# T p*tau_rot p*tau_vib
% # [K] [atm*s] [atm*s]
2000 2.16e-09 2.62e-08
3000 3.02e-09 2.32e-08
4000 3.80e-09 2.19e-08
5000 4.49e-09 2.22e-08
6000 5.41e-09 2.19e-08
8000 6.90e-09 2.19e-08
10000 8.12e-09 2.19e-08
12000 9.42e-09 2.21e-08
15000 1.09e-08 2.15e-08
20000 1.20e-08 2.06e-08
30000 1.39e-08 2.00e-08
40000 1.57e-08 2.03e-08
50000 1.79e-08 2.08e-08
60000 1.99e-08 2.25e-08
80000 2.51e-08 2.45e-08
100000 3.01e-08 2.68e-08    
];

shatalov=[
    0.040679	4.446e-08
0.04735	4.048e-08
0.057686	4.371e-08
0.068977	5.921e-08
0.0792	9.153e-08];
Tshatalov=(2000:500:15000)';
shatalov=[Tshatalov.^(-1/3) ...
    1.5e-12*Tshatalov.^0.5./(1-exp(-2238./Tshatalov))...
    .*exp(86.4*Tshatalov.^(-1/3))];

breen=[
    0.066785	2.9984e-08
0.075672	3.1762e-08
0.099927	2.9525e-08];
Tbreen=(1000:200:3400)';
breen=[Tbreen.^(-1/3) Tbreen./Tbreen*3e-8];

kalogerakis=[
    0.146734	1.3918e-08];

kiefer=[
    0.063121	1.2419e-08
0.069307	2.0487e-08
0.076418	2.6324e-08
0.087304	3.2264e-08];
Tkiefer=(1500:200:3300)';
kiefer=[Tkiefer.^(-1/3) 4.35e-8-7.75e-12*Tkiefer];

ptau_Grover = [%x  ptau_Grover
2975.4  2.0018e-08
4973.2  2.0463e-08
5986.5  2.1335e-08
7998.6  2.2489e-08
9996.5  2.2756e-08
11994.3  2.3026e-08
15019.8  2.6362e-08
];

ptau_Oblapenko = [ % T^(-1/3)  ptau (Oblapenko)
0.030013  1.416e-06
0.03169  1.0664e-06
0.033451  8.116e-07
0.035072  6.583e-07
0.037064  5.149e-07
0.038962  4.1762e-07
0.041092  3.3263e-07
0.043176  2.7476e-07
0.045398  2.249e-07
0.047019  1.9442e-07
0.048824  1.6504e-07
0.051046  1.4138e-07
0.053453  1.2222e-07
0.056184  1.0188e-07
0.058313  9.217e-08
0.060766  8.114e-08
0.063172  7.144e-08
0.065578  6.703e-08
0.068956  5.901e-08
0.072057  5.148e-08
0.075434  4.532e-08
0.078303  4.1754e-08
0.081773  3.7433e-08
0.085891  3.3867e-08
0.090564  2.9277e-08
0.094127  2.6248e-08
0.098152  2.4184e-08
0.102084  2.2281e-08
0.105831  2.0906e-08
0.110088  1.9616e-08
0.115454  1.7747e-08
0.123504  1.5482e-08
0.12961  1.4395e-08
0.1353  1.2788e-08
0.140389  1.1783e-08
0.145015  1.0856e-08
0.149873  1.0002e-08    
];

we=1580.19;
nu=c*we;
theta=h*nu/k;
T=gg(:,3);

M1 = O2;
M2 = O;
kvt10_ssh = zeros(length(T), 1);
kvt10_fho = zeros(length(T), 1);
kvt10_Billing = zeros(length(T), 1);
kvt10_Esposito = zeros(length(T), 1);
for i_T = 1:length(T)
    kvt_temp = kvt_ssh(T(i_T), M1, M2, 1, 1);
    kvt10_ssh(i_T) = kvt_temp(1);
    kvt_temp = kvt_fho_old(T(i_T), M1, M2, 1);
    kvt10_fho(i_T) = kvt_temp(1);
    kvt_temp = kvt_billing(T(i_T), M1, M2, 1, 1);
    kvt10_Billing(i_T) = kvt_temp(1);
    kvt_temp = kvt_Esposito(T(i_T));
    kvt10_Esposito(i_T) = kvt_temp(1);
end
ptau_SSH = 1.363e-22 * T ./ (kvt10_ssh * 1e6 .* (1 - exp(- theta ./ T)));
ptau_FHO = 1.363e-22 * T ./ (kvt10_fho * 1e6 .* (1 - exp(- theta ./ T)));
ptau_Billing = 1.363e-22 * T ./ (kvt10_Billing * 1e6 ...
                                            .* (1 - exp(- theta ./ T)));
ptau_Esposito = 1.363e-22 * T ./ (kvt10_Esposito  ...
                                            .* (1 - exp(- theta ./ T)));

% Millikan−White
mu=1/(1/32+1/16);
mw=10.^(5e-4*mu^0.5*theta^(4/3).*(T.^(-1/3)-0.015*mu^0.25)-8);

% Park correction
sigma_line = 1 * 1e-20; % A^2 -> m2
m_cd = M1.mass * M2.mass / (M1.mass + M2.mass);
ptau_MW_Park_O2O = sqrt(pi * m_cd ./ (4 * k .* T)) / sigma_line ...
    * k .* T * 9.86923e-6;    % s * atm

% Park 1993
a = 47.7;
b = 0.059;
ptau_Park_1993_O2O = exp(a * (T.^(-1/3) - b) - 18.42);

k10FHO=gg(:,5)*1e6;
ptFHO=1.363e-22*T./(k10FHO.*(1-exp(-theta./T)));


figure('Units', 'normalized', 'OuterPosition', [0 0 0.6 0.7]);
semilogy(T.^(-1/3), ptau_Billing, ':','Color', [255,193,7]/255, ...
                                'LineWidth', 2, DisplayName='Billing')
hold on;
semilogy(T.^(-1/3),ptFHO,'-.','Color', [0 0 0],'LineWidth', 2, ...
    DisplayName='FHO RS')
semilogy(T.^(-1/3), ptau_FHO, ':','Color', [0 30 0]/255, ...
    'LineWidth', 2, DisplayName='FHO упрощённые')
semilogy(T.^(-1/3), ptau_SSH, ':','Color', [76,175,80]/255,  ...
    'LineWidth', 3, DisplayName='SSH')
semilogy(T.^(-1/3), ptau_Esposito, 'Color', [0,188,212]/255, ...
                                'LineWidth', 2, DisplayName='Esposito')
semilogy(kiefer(:,1),kiefer(:,2),'o','LineWidth', 1.8, ...
    'color', [33,150,243]/255, 'markerfacecolor',[3,169,244]/255, ...
    DisplayName='Kiefer & Lutz, 1967')
semilogy(breen(:,1),breen(:,2),'x','Color', [0.7 0 0.1],'LineWidth', 2,...
    'markersize', 8, DisplayName='Breen\it et. al.\rm, 1973')
semilogy(kalogerakis(1),kalogerakis(2),'^','LineWidth', 1.5,...
    'color', [211,84,0]/255,'markerfacecolor', [230,126,34]/255,...
    'markersize', 9, DisplayName='Kalogerakis\it et. al.\rm, 2005')
semilogy(shatalov(:,1),shatalov(:,2),'*','Color', [63,81,181]/255,...
    'LineWidth', 1, 'markersize', 7, ...
    DisplayName='Ibraguimova\it et. al.\rm, 2013')
% semilogy(T.^(-1/3), mw, 'LineWidth', 1.5, 'color', [156,39,176]/255, ...
%     DisplayName='Millikan-White')
semilogy(ptau_Torres(:, 1).^(-1/3), ptau_Torres(:, 3), ...
    'Color', [239, 71, 111]/255, 'linewidth', 2, ...
    DisplayName='Torres, 2024')
semilogy(ptau_Grover(:, 1).^(-1/3), ptau_Grover(:, 2), ...
    'sq', 'Color', [17, 138, 178]/255, 'linewidth', 2, ...
    DisplayName='Grover, 2019')
semilogy(ptau_Oblapenko(:, 1), ptau_Oblapenko(:, 2), ...
    'Color', [6, 214, 160]/255, 'linewidth', 2, ...
    DisplayName='Oblapenko, 2018')
% semilogy(T.^(-1/3), ptau_MW_Park_O2O, '--', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='Park поправка')
semilogy(T.^(-1/3), mw + ptau_MW_Park_O2O, '-.', 'LineWidth', 2, ...
    'color', [156,39,176]/255, DisplayName='MW + Park поправка')
% semilogy(T.^(-1/3), ptau_Park_1993_O2O, ':', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='Park 1993')
semilogy(T.^(-1/3), ptau_Park_1993_O2O + ptau_MW_Park_O2O, 'LineWidth', 2.5, ...
    'color', [156,39,176]/255, DisplayName='Park 1993 + Park поправка')
legend('Location', 'east')
% ylim([1e-9 2e-3])
% xlim([0.035 0.153])
% xlim([0.03 0.15])
ylabel('\it{}p\rm\tau, \rmатм \cdot c');
xlabel('\it{}T\rm ^{-1/3}, K');
grid on;
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% set(gcf, 'color', 'none');
%  set(gca, 'color', 'none');

%% 2) p*tau for NO-NO
we = NO.we(1)/100;
nu = c * we;
theta = h * nu / k;
T = gg(:, 3);

M1 = NO;
M2 = NO;
kvt10_ssh = zeros(length(T), 1);
kvt10_fho = zeros(length(T), 1);
for i_T = 1:length(T)
    kvt_temp = kvt_ssh(T(i_T), M1, M2, 1, 1);
    kvt10_ssh(i_T) = kvt_temp(1);
    kvt_temp = kvt_fho_old(T(i_T), M1, M2, 1);
    kvt10_fho(i_T) = kvt_temp(1);
end
ptau_SSH = 1.363e-22 * T ./ (kvt10_ssh * 1e6 .* (1 - exp(- theta ./ T)));
ptau_FHO = 1.363e-22 * T ./ (kvt10_fho * 1e6 .* (1 - exp(- theta ./ T)));

ptau_Torres_NONO = [ % T^-1/3	ptau_Torres_NONO
0.0215946	2.0043e-07
0.0231995	1.8574e-07
0.0255894	1.791e-07
0.0271945	1.7443e-07
0.029317	1.7558e-07
0.0322601	1.8029e-07
0.0368979	2.0176e-07
0.0405546	2.1988e-07
0.0437119	2.4282e-07
0.0464947	2.6994e-07
0.0500806	3.336e-07
0.0550406	5.0614e-07
0.0585914	7.3072e-07
0.0630349	1.3432e-06
];

ptau_Oblapenko_NONO = [ % T^-1/3  ptau_Oblapenko_NO
0.0299826  1.0629e-06
0.0313449  9.386e-07
0.0328886  8.417e-07
0.0348259  7.205e-07
0.0370959  6.262e-07
0.0393961  5.528e-07
0.0415749  5.113e-07
0.0437539  4.62e-07
0.0461445  4.272e-07
0.0490494  3.9192e-07
0.052529  3.7081e-07
0.0567951  3.5075e-07
0.0645701  3.5532e-07
0.072678  3.5713e-07
0.0793638  3.6757e-07
0.0851724  3.7259e-07
0.0933106  3.716e-07
0.099089  3.709e-07
0.1062293  3.4511e-07
0.1125825  3.4439e-07
0.1205093  3.279e-07
0.1280427  3.1467e-07
0.1345473  3.0917e-07
0.140084  2.9689e-07
0.1457416  2.8729e-07
0.1501588  2.8029e-07
];

ptau_Moser_NONO = [ % T^-1/3	ptau Moser (Oblapenko)
0.0659486	2.2096e-06
0.0737546	2.0239e-06
0.0794428	1.8551e-06
];

ptau_Glanzer_1975_NONO = [ % T^-1/3  Glanzer_1975 (Oblapenko)
0.071926  2.767e-07
0.0748025  3.2952e-07
0.0781564  3.4901e-07
0.082332  2.7175e-07
0.087458  3.0056e-07
0.094224  2.8983e-07
0.1036952  3.2641e-07
];

ptau_Glanzer_1977_NONO = [ % T^-1/3  Glanzer_1977 (Oblapenko)
0.0694208  7.123e-07
0.0705558  6.641e-07
0.0717663  6.288e-07
0.0731131  5.908e-07
0.0745203  5.55e-07
0.0760182  5.255e-07
0.0776373  4.898e-07
0.0794378  4.601e-07
0.081435  4.322e-07
0.0836138  4.0281e-07
0.0860646  3.8719e-07
0.0887728  3.608e-07
0.0919043  3.454e-07
0.0955047  3.3576e-07
0.0997856  3.3014e-07
0.1049741  3.3214e-07
0.1114329  3.4586e-07
0.1198577  3.8293e-07
0.1315491  4.72e-07
0.1498062  5.352e-07
];

ptau_Hancock_NONO = [ % T^-1/3  Hancock (Oblapenko)
0.1502931  4.47e-08
];

ptau_Kamimoto_NONO_Oblapenko = [ % T^-1/3   Kamimoto (Oblapenko, wrong)
0.0694388  4.762e-07
0.07136  4.722e-07
0.0734776  4.792e-07
0.0758222  4.845e-07
0.0785901  5.012e-07
0.0817968  5.144e-07
0.0855934  5.259e-07
0.0902221  5.353e-07
0.0960305  5.49e-07
0.1036237  5.826e-07
];

ptau_Horiguchi_NONO = [ % T^-1/3  Horiguchi (Oblapenko)
0.1503287  3.349e-08
];

ptau_Wray_NONO = [ % T^-1/3  Wray (Oblapenko)
0.0501815  2.2215e-09
0.0507812  2.3985e-09
0.051441  2.5895e-09
0.0521247  2.8129e-09
0.0528685  3.0556e-09
0.0536363  3.2989e-09
0.0544519  3.6952e-09
0.0553278  4.0262e-09
0.0562393  4.566e-09
0.0572469  5.161e-09
0.0582905  5.907e-09
0.0594662  6.697e-09
0.0607497  7.807e-09
0.0620932  9.156e-09
0.0635808  1.0771e-08
0.0652241  1.3349e-08
0.0670836  1.6292e-08
0.0691591  2.0189e-08
0.0715225  2.5797e-08
0.074222  3.3368e-08
0.077426  4.198e-08
0.0812543  5.215e-08
0.0859233  6.339e-08
0.0919493  7.287e-08
0.1001125  8.346e-08
];

T_Streicher = [1900; T(10:23); 8600];
ptau_Streicher = 1.1e-3 * exp(70 * T_Streicher.^(-1/3)) / 1e6;

T_Kamimoto = T(7:12);
ptau_KAMIMOTO = 10.^(6.1 * T_Kamimoto .^ (-1/3) - 1.1) / 1e6;

T_Breshears = T(3:10);
ptau_Breshears = 3.48 * ...
    exp(-38.6 * T_Breshears.^(-1/3) + 1582 * T_Breshears .^(-1)) / 1e6;

% Millikan−White
mu=1/(1/30.01+1/30.01);
mw_NONO=10.^(5e-4*mu^0.5*theta^(4/3).*(T.^(-1/3)-0.015*mu^0.25)-8);

% Park correction
sigma_line = 1 * 1e-20; % A^2 -> m2
m_cd = M1.mass * M2.mass / (M1.mass + M2.mass);
ptau_MW_Park_NONO = sqrt(pi * m_cd ./ (4 * k .* T)) / sigma_line ...
    * k .* T * 9.86923e-6;    % s * atm

% Park 1993
a = 49.5;
b = 0.042;
ptau_Park_1993_NONO = exp(a * (T.^(-1/3) - b) - 18.42);


figure('Units', 'normalized', 'OuterPosition', [0 0 0.6 0.7]);
semilogy(T.^(-1/3), ptau_FHO, ':','Color', [0 30 0]/255, ...
    'LineWidth', 1.5, DisplayName='FHO')
hold on;
semilogy(T.^(-1/3), ptau_SSH, ':','Color', [76,175,80]/255,  ...
    'LineWidth', 3, DisplayName='SSH')
semilogy(T_Streicher.^(-1/3), ptau_Streicher, '-.', ...
    'Color', [255, 209, 102]/255, 'LineWidth', 2, DisplayName='Streicher')
semilogy(T_Kamimoto.^(-1/3), ptau_KAMIMOTO, '--', 'LineWidth', 3, ...
    'Color', [17, 138, 178]/255, DisplayName='Kamimoto, 1970')
semilogy(ptau_Oblapenko_NONO(:, 1), ptau_Oblapenko_NONO(:, 2), ...
    'Color', [6, 214, 160]/255, 'LineWidth', 1.5, ...
    DisplayName='Oblapenko 2018')
semilogy(ptau_Torres_NONO(:, 1), ptau_Torres_NONO(:, 2), ...
    'Color', [239, 71, 111]/255, 'LineWidth', 2, ...
    DisplayName='Torres, 2024')
semilogy(T_Breshears.^(-1/3), ptau_Breshears, 'LineWidth', 2, ...
    'Color', [142, 68, 173]/255, DisplayName='Breshears, 1969')
semilogy(ptau_Moser_NONO(:,1), ptau_Moser_NONO(:,2), 'o', ...
    'LineWidth', 1.8, 'color', [251, 111, 146]/255, ...
    'markerfacecolor', [251, 111, 146]/255, DisplayName='Moser')
semilogy(ptau_Glanzer_1975_NONO(:, 1), ptau_Glanzer_1975_NONO(:, 2), ...
    '^', 'LineWidth', 1.5, 'color', [211,84,0]/255, ...
    'markerfacecolor', [230,126,34]/255, 'markersize', 6, ...
    DisplayName='Glanzer, 1975')
semilogy(ptau_Glanzer_1977_NONO(:, 1), ptau_Glanzer_1977_NONO(:, 2), ...
    'v', 'LineWidth', 1.5, 'color', [211,84,0]/255, ...
    'markerfacecolor', [230,126,34]/255, 'markersize', 6, ...
    DisplayName='Glanzer, 1977')
semilogy(ptau_Hancock_NONO(:, 1), ptau_Hancock_NONO(:, 2), 'x', ...
    'Color', [0.7 0 0.1],'LineWidth', 2, 'markersize', 8, ...
    DisplayName='Hancock')
semilogy(ptau_Horiguchi_NONO(:,1), ptau_Horiguchi_NONO(:,2), '*', ...
    'Color', [63,81,181]/255, 'LineWidth', 1, 'markersize', 7, ...
    DisplayName='Horiguchi')
% semilogy(T.^(-1/3), mw_NONO, 'LineWidth', 1.5, ...
%     'color', [156,39,176]/255, DisplayName='Millikan-White')
semilogy(ptau_Wray_NONO(:, 1), ptau_Wray_NONO(:, 2), ...
    'sq', 'Color', [17, 138, 178]/255, 'linewidth', 2, ...
    DisplayName='Wray, 1962')
% semilogy(T.^(-1/3), ptau_MW_Park_NONO, '--', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='Park поправка')
semilogy(T.^(-1/3), mw_NONO + ptau_MW_Park_NONO, '-.', 'LineWidth', 2, ...
    'color', [156,39,176]/255, DisplayName='MW + Park поправка')
% semilogy(T.^(-1/3), ptau_Park_1993_NONO, ':', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='Park 1993')
semilogy(T.^(-1/3), ptau_Park_1993_NONO + ptau_MW_Park_NONO, 'LineWidth', 2.5, ...
    'color', [156,39,176]/255, DisplayName='Park 1993 + Park поправка')
% legend('FHO', 'SSH', 'Streicher', 'Kamimoto, 1970', 'Oblapenko 2018', ...
%     'Torres, 2024', 'Breshears, 1969', 'Moser', 'Glanzer, 1975', ...
%     'Glanzer, 1977', 'Hancock', 'Horiguchi', ...'Kamimoto Oblapenko', ...
%     'Millikan-White', 'Wray', ...
%     'Location', 'east')
legend('Location', 'east')
ylim([1e-9 1e1])
% xlim([0.035 0.153])
% xlim([0.03 0.15])
ylabel('\it{}p\rm\tau, \rmатм \cdot c');
xlabel('\it{}T\rm ^{-1/3}, K');
grid on;
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% set(gcf, 'color', 'none');
%  set(gca, 'color', 'none');

%% 3) p*tau for CO-CO
we = CO.we(1)/100;
nu = c * we;
theta = h * nu / k;
T = gg(:, 3);

M1 = CO;
M2 = CO;
kvt10_ssh = zeros(length(T), 1);
kvt10_fho = zeros(length(T), 1);
kvt10_DeLeon = zeros(length(T), 1);
ptau_DeLeon_clear = zeros(length(T), 1);
for i_T = 1:length(T)
    kvt_temp = kvt_ssh(T(i_T), M1, M2, 1, 1);
    kvt10_ssh(i_T) = kvt_temp(1);
    kvt_temp = kvt_fho_old(T(i_T), M1, M2, 1);
    kvt10_fho(i_T) = kvt_temp(1);
    kvt_temp = kvt_DeLeon_CO(T(i_T), M1);
    kvt10_DeLeon(i_T) = kvt_temp(1);
    ptau_DeLeon_clear(i_T) = ptau_DeLeon_CO(T(i_T), M1) / 1e6;
end
ptau_SSH_COCO = 1.363e-22 * T ./ ...
    (kvt10_ssh * 1e6 .* (1 - exp(- theta ./ T)));
ptau_FHO_COCO = 1.363e-22 * T ./ ...
    (kvt10_fho * 1e6 .* (1 - exp(- theta ./ T)));
ptau_DeLeon_COCO = 1.363e-22 * T ./ ...
    (kvt10_DeLeon * 1e6 .* (1 - exp(- theta ./ T)));

Hooker=[
1367.33	0.00030435
1425.26	0.00024758
1512.68	0.00020729
1608.63	0.0002062
1653	0.00017404
1676.17	9.887e-05
1722.45	0.000133967
1854.69	0.00010608
1968.04	9.9475e-05
2049.17	7.5163e-05
2101.39	6.3901e-05
2205.29	4.0174e-05
2371.8	3.0424e-05
2372.3	3.4466e-05
2633.89	2.2751e-05
2647.06	1.9931e-05
];
Matthews=[
2196.96	5.778e-05
2293.85	2.7489e-05
2414.03	3.2717e-05
2791.1	1.59536e-05
];
Gaydon=[
2184.1	4.3889e-05
2336.46	2.4426e-05
2422.24	2.209e-05
2544.09	2.452e-05
2632.94	2.9199e-05
];
Russo=[
1432.41	0.0002091
1602.18	0.000110108
1842.69	6.1414e-05
1942.64	3.8332e-05
1958.56	4.0788e-05
2231.65	3.0946e-05
2425.71	1.9786e-05
2534.91	1.7562e-05
];

deleon=[
300  2.29395e-27
400  1.66563e-26
500  7.93763e-26
600  2.85992e-25
700  8.45683e-25
800  2.15985e-24
900  4.92759e-24
1000  1.02796e-23
1500  1.68518e-22
2000  1.17309e-21
2500  5.11471e-21
3000  1.66278e-20
3500  4.42349e-20
4000  1.01777e-19
4500  2.09846e-19
5000  3.97179e-19
5500  7.01963e-19
6000  1.17305e-18
6500  1.87107e-18
7000  2.86943e-18
7500  4.25519e-18
8000  6.12989e-18
8500  8.61021e-18
9000  1.18286e-17
9500  1.59338e-17
10000  2.10914e-17
11000  3.53112e-17
12000  5.61595e-17
13000  8.55938e-17
14000  1.25862e-16
15000  1.79505e-16
16000  2.49355e-16
17000  3.38537e-16
18000  4.50463e-16
19000  5.8883e-16
20000  7.57616e-16];
deleon(:,2)=deleon(:,2)*1e6;
ptDL=1.363e-22*deleon(:,1)./(deleon(:,2).*(1-exp(-theta./deleon(:,1))));


% Millikan−White
mu=1/(1/28.01+1/28.01);
mw_COCO = 10.^(5e-4*mu^0.5*theta^(4/3).*(T.^(-1/3)-0.015*mu^0.25)-8);

% Park correction
sigma_line = 1 * 1e-20; % A^2 -> m2
m_cd = M1.mass * M2.mass / (M1.mass + M2.mass);
ptau_MW_Park_COCO = sqrt(pi * m_cd ./ (4 * k .* T)) / sigma_line ...
    * k .* T * 9.86923e-6;    % s * atm

% Park 1993
a = 198;
b = 0.0290;
ptau_Park_1994_COCO = exp(a * (T.^(-1/3) - b) - 18.42);

figure('Units', 'normalized', 'OuterPosition', [0 0 0.6 0.7]);
semilogy(T.^(-1/3), ptau_FHO_COCO, ':','Color', [0 30 0]/255, ...
    'LineWidth', 1.5, DisplayName='FHO')
hold on;
semilogy(T.^(-1/3), ptau_SSH_COCO, ':','Color', [76,175,80]/255,  ...
    'LineWidth', 3, DisplayName='SSH')
semilogy(T.^(-1/3), ptau_DeLeon_COCO, ':','Color', [227, 178, 60]/255,  ...
    'LineWidth', 3, DisplayName='DeLeon')
% semilogy(deleon(:,1).^(-1/3), ptDL, 'Color', [227, 178, 60]/255, ...
%     DisplayName='DeLeon')
% semilogy(T.^(-1/3), ptau_DeLeon_clear, '--','Color', [227, 178, 60]/255,  ...
%     'LineWidth', 3, DisplayName='DeLeon ptau')
% hold on;
semilogy(Hooker(:,1).^(-1/3), Hooker(:,2), '*', DisplayName='Hooker, exp')
semilogy(Matthews(:,1).^(-1/3), Matthews(:,2), '^', ...
    DisplayName='Matthews, exp')
semilogy(Gaydon(:,1).^(-1/3), Gaydon(:,2), 'o', ...
    DisplayName='Gaydon, exp')
semilogy(Russo(:,1).^(-1/3), Russo(:,2), 'sq', DisplayName='Russo, exp')
% semilogy(T.^(-1/3), mw_COCO, 'LineWidth', 1.5, ...
%     'color', [156,39,176]/255, DisplayName='Millikan-White')
% semilogy(T.^(-1/3), mw_COCO + ptau_MW_Park_COCO, '-.', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='MW + Park поправка')
% semilogy(T.^(-1/3), ptau_Park_1994_COCO, ':', 'LineWidth', 2, ...
%     'color', [156,39,176]/255, DisplayName='Park 1994')
semilogy(T.^(-1/3), ptau_Park_1994_COCO + ptau_MW_Park_COCO, 'LineWidth', 2.5, ...
    'color', [156,39,176]/255, DisplayName='Park 1994 + Park поправка')
 % legend('FHO', 'SSH', 'DeLeon', 'Hooker, exp', 'Matthews, exp', ...
 %     'Gaydon, exp', 'Russo, exp', 'Location','SouthEast')
legend('Location','SouthEast')
ylabel('p\tau, sec. atm')
xlabel('T^{-1/3}, K')
xlim([0.02 0.16])
ylim([0.99e-7 1e1]);
grid on;
title('CO + CO');
%%

rmpath('../src/')