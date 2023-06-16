%addpath('..\data\Streicher experiment\');
addpath('..\data\NO experiment\');
load('NO_behind_ReflSW_withexch.mat');
load('NO_between_SWs_withexch.mat');
info=["2% No.1", "2% No.2", "2% No.3" ,"2% No.4" ,"1% No.5", "1% No.6",...
    "1% No.7","1% No.8","0.4% No.9","0.4% No.10","0.4% No.11"];
for var=[1:11]
%var=6; %1 - 2-02 T=3560 P=0.561;  2 - 2-14 T=5460 P=0.325; 3 - 2-32 T=7070
%P=0.119; 4 - 2-38 T=8730 P=0.137
% 5 - 1-01 T=2360 P=0.757; 6 - 1-04 T=3470 P=0.584; 7 - 1-08 T=5580K P=0.252;
% 8 - 1-12 T=7090K P=0.274; 9 - 04-02 T=2180K P=0.952; 10 - 04-11 T=3470K 
%P=0.668;  11 - 04-22 T=5650K P=0.545;
i_vibr=1;
rel=2;
switch var
    case 1
    TT=readmatrix('T_2%_1.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_2%_1.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=600;
    case 2
    TT=readmatrix('T_2%_2.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_2%_2.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=300;
    case 3
    TT=readmatrix('T_2%_3.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_2%_3.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=150;
    case 4
    TT=readmatrix('T_2%_4.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_2%_4.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=70;
    case 5
    TT=readmatrix('T_1%_1.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_1%_1.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=600;
    case 6
    TT=readmatrix('T_1%_2.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,3);
    n_exp=readmatrix('n_1%_2.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=600;
    case 7
    Test=readmatrix("T_1%_3.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,2);
    Tv_exp=Test(:,3);
    n_exp=readmatrix('n_1%_3.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=300;
    case 8
    Test=readmatrix("T_1%_4.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,2);
    Tv_exp=Test(:,3);
    n_exp=readmatrix('n_1%_3.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=150;
    case 9
    Test=readmatrix("T_0_4%_1.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,2);
    Tv_exp=Test(:,3);
    n_exp=readmatrix('n_0_4%_1.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=600;
    case 10
    Test=readmatrix("T_0_4%_2.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,2);
    Tv_exp=Test(:,3);
    n_exp=readmatrix('n_0_4%_2.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=600;
    case 11
    Test=readmatrix("T_0_4%_3.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,2);
    Tv_exp=Test(:,3);    
    n_exp=readmatrix('n_0_4%_3.csv');
    time_n_exp=n_exp(:,1);
    n_exp=n_exp(:,2);
    tlim=300;
end

figure("Position", [0, 0, 900, 800])
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "tight")
nexttile
hold on
title("Case " + info(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
plot(time_T_exp, T_exp, 'k-', 'LineWidth', 2, 'DisplayName', "T - experiment ");
legend('Location','se');
xlim([-10 300]);
ylim([0 max(T_exp)+500]);
xlabel("t, \mu s");
ylabel("T, K");
hold off
grid minor
nexttile
hold on
title("Case " + info(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
plot(time_T_exp, Tv_exp, 'k-', 'LineWidth', 2, 'DisplayName', "Tv - experiment");
xlim([-10 300]);
%ylim([0 6000]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("T_v, K");
hold off
grid minor
nexttile
hold on
title("Case " + info(var));
plot(time_n_exp, n_exp, 'k-', 'LineWidth', 2, 'DisplayName', "n_{NO} - experiment" + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).nNO*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n_{NO} - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).nNO*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n_{NO} - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).nNO*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n_{NO} - U=inf " );
xlim([-10 tlim]);
ylim([0 max(dat1(i_vibr,2,var,rel).nNO*1e3)+4]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("n_{NO}, mmol/m^3");
hold off
grid minor
nexttile
hold on
title("Case " + info(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
xlim([-10 10000]);
%ylim([0 6000]);
legend('Location','se');
xlabel("t, \mu s");
ylabel("T, K");
hold off
grid minor
end;