%%Comparison with experimental data and plots
%Streicher's 2021 experiment with O2\O\Ar mixtures
clearvars;
load('..\data\O2_Ar Streicher experiment\output_ReflSW_O2_Ar.mat');
load('..\data\O2_Ar Streicher experiment\O2_Ar_Streicher21_experiment2.mat');
info=["50% No.1 (03)", "50% No.2 (11)", "50% No.3 (14)" ,"20% No.1 (02)" ,"20% No.2 (08)", "20% No.3 (14)",...
    "100% No.1 (01)","100% No.2 (06)","100% No.3 (08)"];

%Calculated data
dat=data_betweenSWs_O2Ar;
dat1=data_behindRSW_O2Ar;

for var=1:9
%testcases
%var: %1 - 50-03 T=8110 P=75;  2 - 50-11 T=10470 P=53; 3 - 50-13 T=11410 P=30; 4 - 20-02 T=7840 P=130
% 5 - 20-08 T=10310 P=97; 6 - 20-14 T=13830 P=33; 7 - 100-01 T=6230K P=57;
% 8 - 100-06 T=7940K P=41; 9 - 100-08 T=9560K P=34;
i_vibr=1; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on

%initialization of experimental data for temperature
% time_T_exp=data_experiment(var).T(:,1);
% T_exp=data_experiment(var).T(:,2);
% Tv_exp=data_experiment(var).T(:,3);
% %initialization of experimental data for number density
% time_n_exp=data_experiment(var).n(:,1);
% nm_exp=data_experiment(var).n(:,2);


time_Tv_err=[];
Tv_err=[];
err_Tv=[];
time_n_err=[];
n_err=[];
err_n=[];
%initialization of experimental data for temperature
j=1;
while (data_experiment(var).Tv(j,3)~=0)
j=j+1;
end
time_Tv_exp=data_experiment(var).Tv(1:j-1,1);
Tv_exp=data_experiment(var).Tv(1:j-1,2);
[time_Tv_exp, I]=sort(time_Tv_exp);
Tv_exp=Tv_exp(I);
tlim_T=max(time_Tv_exp);
for i=j+1:length(data_experiment(var).Tv(:,3))-1
    if (data_experiment(var).Tv(i,3)~=-1)&&(data_experiment(var).Tv(i-1,3)==data_experiment(var).Tv(i,3))&&(data_experiment(var).Tv(i+1,3)==data_experiment(var).Tv(i,3))
        time_Tv_err=[time_Tv_err, (data_experiment(var).Tv(i,1)+data_experiment(var).Tv(i-1,1)+data_experiment(var).Tv(i+1,1))/3];
        Tv_err=[Tv_err, data_experiment(var).Tv(i-1,2)];
        err_Tv=[err_Tv, (data_experiment(var).Tv(i,2) - data_experiment(var).Tv(i+1,2))/2];
    end
end

time_T_exp=data_experiment(var).T(:,1);
T_exp=data_experiment(var).T(:,2);
[time_T_exp, I]=sort(time_T_exp);
T_exp=T_exp(I);

j=1;
while (data_experiment(var).n(j,3)~=0)
j=j+1;
end
time_n_exp=data_experiment(var).n(1:j-1,1);
tlim_n=max(time_n_exp);
n_exp=data_experiment(var).n(1:j-1,2);
[time_n_exp, I]=sort(time_n_exp);
n_exp=n_exp(I);
for i=(j+1):(length(data_experiment(var).n(:,3))-1)
    if (data_experiment(var).n(i,3)~=-1)&&(data_experiment(var).n(i-1,3)==data_experiment(var).n(i,3))&&(data_experiment(var).n(i+1,3)==data_experiment(var).n(i,3))
        time_n_err=[time_n_err, (data_experiment(var).n(i,1)+data_experiment(var).n(i-1,1)+data_experiment(var).n(i+1,1))/3];
        n_err=[n_err, data_experiment(var).n(i-1,2)];
        err_n=[err_n, (data_experiment(var).n(i,2) - data_experiment(var).n(i+1,2))/2];
    end
end



% if var==1
% time_p_exp=data_experiment(var).p(:,1);
% p_exp=data_experiment(var).p(:,2);
% end

%%Temperature & Number density
figure("Position", [0, 0, 900, 800])
t=tiledlayout(2, 2, "TileSpacing", "compact");
title(t, "Case " + info(var));
nexttile
hold on
plot(time_T_exp, T_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "T - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
legend('Location','best');
xlim([0 50]);
ylim([0 14000]);
xlabel("t, \mus");
ylabel("T, K");
hold off
grid minor
nexttile
hold on
title("Case " + info(var));
plot(time_Tv_exp, Tv_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', "Tv - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
errorbar(time_Tv_err, Tv_err, err_Tv, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
xlim([0 50]);
legend('Location','best');
xlabel("t, \mus");
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
ylim([0 11000]);
legend('Location','best');
xlabel("t, \mus");
ylabel("T_v, K");
hold off
grid minor

nexttile 
hold on
plot(time_n_exp, n_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', "n_m - exp case " + num2str(var));
plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).nm_n*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=D/6k " );
plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).nm_n*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n_m - U=3T " );
plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).nm_n*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n_m - U=inf " );
errorbar(time_n_err, n_err, err_n, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
legend('Location','best');
xlabel("t, \mus");
ylabel("n_{O_2}, mmol/m^3");
xlim([0 50]);
ylim([0 max(dat1(i_vibr,2,var,rel).nm_n*1e3) + 4]);
hold off
grid minor
end

%%Pressure
% if var==1
% PPP(:,1)=[75 53 30 130 97 33 57 41 34];
% PPP(:,2)=[0.12 0.30 0.36 0.09 0.23 0.29 0.10 0.15 0.24];
% figure
% hold on
% plot(time_p_exp, p_exp, 'DisplayName',"p - raw data");
% plot(time_p_exp(16:end), time_p_exp(16:end)*PPP(1,2)+PPP(1,1), 'k-','LineWidth', 1.5, 'DisplayName',"p - interpolated data");
% plot(dat1(i_vibr,2,1,rel).time, dat1(i_vibr,2,1,rel).p,'r-', 'LineWidth', 1.5, 'DisplayName', "p - U=D/6k SSH");
% plot(dat1(i_vibr,3,1,rel).time, dat1(i_vibr,3,1,rel).p,'b-', 'LineWidth', 1.5, 'DisplayName', "p - U=3T SSH");
% plot(dat1(i_vibr,4,1,rel).time, dat1(i_vibr,4,1,rel).p,'m-', 'LineWidth', 1.5, 'DisplayName', "p - U=\infty SSH");
% legend('Location','best');
% xlabel("t, \mus");
% ylabel("p, Torr");
% xlim([-20 100]);
% ylim([-10 100]);
% hold off;
% grid minor;
% end
% 
