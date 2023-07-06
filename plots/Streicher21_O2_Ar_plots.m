%%Сравнение с экспериментом, графики
addpath('..\data\O2_Ar Streicher experiment\');
load('O2Ar_Streicher_behind_ReflSW.mat');
load('O2Ar_Streicher_between_SWs.mat');
load('O2_Ar_Streicher21_experiment.mat');
info=["50% No.1 (03)", "50% No.2 (11)", "50% No.3 (14)" ,"20% No.1 (02)" ,"20% No.2 (08)", "20% No.3 (14)",...
    "100% No.1 (01)","100% No.2 (06)","100% No.3 (08)"];
for var=[1:9]
%testcases
%var: %1 - 50-03 T=8110 P=75;  2 - 50-11 T=10470 P=53; 3 - 50-13 T=11410 P=30; 4 - 20-02 T=7840 P=130
% 5 - 20-08 T=10310 P=97; 6 - 20-14 T=13830 P=33; 7 - 100-01 T=6230K P=57;
% 8 - 100-06 T=7940K P=41; 9 - 100-08 T=9560K P=34;
i_vibr=1; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=1; %switcher of relaxation between SWs: 1 - off, 2 - on

%initialization of experimental data for temperature
time_T_exp=data_experiment(var).T(:,1);
T_exp=data_experiment(var).T(:,2);
Tv_exp=data_experiment(var).T(:,3);
%initialization of experimental data for number density
time_n_exp=data_experiment(var).n(:,1);
nm_exp=data_experiment(var).n(:,2);

if var==1
time_p_exp=readmatrix('p.csv');
p_exp=time_p_exp(:,2);
time_p_exp=time_p_exp(:,1);
end

%%Temperature & Number density
figure("Position", [0, 0, 900, 800])
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "tight")
nexttile
hold on
plot(time_T_exp, T_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "T - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
legend('Location','ne');
xlim([0 50]);
xlabel("t, \mu s");
ylabel("T, K");
hold off
grid minor
nexttile
hold on
title("Case " + info(var));
plot(time_T_exp, Tv_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', "Tv - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
xlim([0 50]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("T_v, K");
hold off
grid minor

nexttile
hold on
title("Case " + info(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
xlim([0 450]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("T_v, K");
hold off
grid minor

nexttile 
hold on
plot(time_n_exp, nm_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "n_m - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).nm_n*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).nm_n*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).nm_n*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n_m - U=inf " );
legend('Location','ne');
xlabel("t, \mu s");
ylabel("n_{O_2}, mmol/m^3");
xlim([0 50]);
hold off
grid minor


%%Pressure

% PPP(:,1)=[75 53 30 130 97 33 57 41 34];
% PPP(:,2)=[0.12 0.30 0.36 0.09 0.23 0.29 0.10 0.15 0.24];
% figure
% hold on
% plot(time_p_exp, p_exp, 'DisplayName',"p - raw data [9]");
% plot(time_p_exp(16:end), time_p_exp(16:end)*PPP(1,2)+PPP(1,1), 'k-','LineWidth', 1.5, 'DisplayName',"p - interpolated data [9]");
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).p,'r-', 'LineWidth', 1.5, 'DisplayName', "p - U=D/6k SSH");
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).p,'b-', 'LineWidth', 1.5, 'DisplayName', "p - U=3T SSH");
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).p,'m-', 'LineWidth', 1.5, 'DisplayName', "p - U=\infty SSH");
% legend('Location','se');
% xlabel("t, \mus");
% ylabel("p, Torr");
% xlim([-20 100]);
% ylim([-10 100]);
% hold off;
% grid minor;
% 

% PPP(var, 3)=dat1(i_vibr,2,var,rel).p(1);
% a=polyfit(dat1(i_vibr,2,var,rel).time(1:118),dat1(i_vibr,2,var,rel).p(1:118), 1);
% PPP(var, 4)=a(1);
% PPP(var, 5)=a(2);
% Pressure=array2table(PPP, "VariableNames",["p0 - exp","dp/dt - exp", "p_behindRSW U=D/6k",...
%     "dp/dt U=D/6k", "p0_interp U=D/6k"]);
end


