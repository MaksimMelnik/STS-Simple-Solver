%%Сравнение с экспериментом, графики
addpath('..\data\Streicher experiment\');
load('O2Ar_Streicher_behind_ReflSW.mat');
info=["50% No.1 (03)", "50% No.2 (11)", "50% No.3 (14)" ,"20% No.1 (02)" ,"20% No.2 (08)", "20% No.3 (14)",];
for var=[1:6]
%var=6; %1 - 50-03 T=8110 P=75;  2 - 50-11 T=10470 P=53; 3 - 50-13 T=11410 P=30; 4 - 20-02 T=7840 P=130
% 5 - 20-08 T=10310 P=97; 6 - 20-14 T=13830 P=33
i_vibr=1;

switch var
    case 1
    TT=readmatrix('T_50.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,5);
    nm_exp=readmatrix('nm_50_1.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    case 2
    TT=readmatrix('T_50.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,3);
    Tv_exp=TT(:,6);
    nm_exp=readmatrix('nm_50_2.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    case 3
    TT=readmatrix('T_50.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,4);
    Tv_exp=TT(:,7);
    nm_exp=readmatrix('nm_50_3.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    case 4
    TT=readmatrix('T_20.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,5);
    nm_exp=readmatrix('nm_20_1.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    case 5
    TT=readmatrix('T_20.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,3);
    Tv_exp=TT(:,6);
    nm_exp=readmatrix('nm_20_2.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    case 6
    TT=readmatrix('T_20.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,4);
    Tv_exp=TT(:,7);
    nm_exp=readmatrix('nm_20_3.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
end

%%Temperature
figure("Position", [0, 0, 1300, 400])
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "tight")
nexttile
hold on
plot(time_T_exp, T_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "T - exp case " + num2str(var));
plot(dat1(i_vibr,2,var).time, dat1(i_vibr,2,var).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var).time, dat1(i_vibr,3,var).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var).time, dat1(i_vibr,4,var).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
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
plot(dat1(i_vibr,2,var).time, dat1(i_vibr,2,var).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var).time, dat1(i_vibr,3,var).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var).time, dat1(i_vibr,4,var).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
xlim([0 50]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("T_v, K");
hold off
grid minor
nexttile 
hold on
plot(time_n_exp, nm_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "n_m - exp case " + num2str(var));
plot(dat1(i_vibr,2,var).time, dat1(i_vibr,2,var).nm_n*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=D/6k " );
plot(dat1(i_vibr,3,var).time, dat1(i_vibr,3,var).nm_n*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=3T " );
plot(dat1(i_vibr,4,var).time, dat1(i_vibr,4,var).nm_n*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n_m - U=inf " );
legend('Location','ne');
xlabel("t, \mu s");
ylabel("n_{O_2}, mmol/m^3");
xlim([0 50]);
hold off
grid minor


end