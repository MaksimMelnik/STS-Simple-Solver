%%Comparison with experimental data and plots
%Streicher's 2022 experiment with NO\O\N\O2\N2\Ar mixtures

load('..\data\NO Streicher experiment\NO_behind_ReflSW_withexch.mat');
load('..\data\NO Streicher experiment\NO_between_SWs_withexch.mat');
load('..\data\NO Streicher experiment\NO_Streicher22_experiment.mat');
info=["2% No.1", "2% No.2", "2% No.3" ,"2% No.4" ,"1% No.5", "1% No.6",...
    "1% No.7","1% No.8","0.4% No.9","0.4% No.10","0.4% No.11"];
for var=1:11
%testcases     
%var: %1 - 2-02 T=3560 P=0.561;  2 - 2-14 T=5460 P=0.325; 3 - 2-32 T=7070
%P=0.119; 4 - 2-38 T=8730 P=0.137
% 5 - 1-01 T=2360 P=0.757; 6 - 1-04 T=3470 P=0.584; 7 - 1-08 T=5580K P=0.252;
% 8 - 1-12 T=7090K P=0.274; 9 - 04-02 T=2180K P=0.952; 10 - 04-11 T=3470K 
%P=0.668;  11 - 04-22 T=5650K P=0.545;
i_vibr=2; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on

%initialization of experimental data for temperature
time_T_exp=data_experiment(var).T(:,1);
T_exp=data_experiment(var).T(:,2);
Tv_exp=data_experiment(var).T(:,3);

%initialization of experimental data for number density
time_n_exp=data_experiment(var).n(:,1);
n_exp=data_experiment(var).n(:,2);
tlim=max(time_n_exp);

err=0.045;
   

figure("Position", [0, 0, 900, 800])
tiledlayout(2, 2, "TileSpacing", "compact")
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
errorbar(time_n_exp, n_exp, n_exp*err, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1)
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
end