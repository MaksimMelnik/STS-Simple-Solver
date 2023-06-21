%%Сравнение с экспериментом, графики
addpath('..\data\Streicher experiment\');
load('O2Ar_Streicher_behind_ReflSW.mat');
load('O2Ar_Streicher_between_SWs.mat');
info=["50% No.1 (03)", "50% No.2 (11)", "50% No.3 (14)" ,"20% No.1 (02)" ,"20% No.2 (08)", "20% No.3 (14)",...
    "100% No.1 (01)","100% No.2 (06)","100% No.3 (08)"];
for var=[1:9]
%var=6; %1 - 50-03 T=8110 P=75;  2 - 50-11 T=10470 P=53; 3 - 50-13 T=11410 P=30; 4 - 20-02 T=7840 P=130
% 5 - 20-08 T=10310 P=97; 6 - 20-14 T=13830 P=33; 7 - 100-01 T=6230K P=57;
% 8 - 100-06 T=7940K P=41; 9 - 100-08 T=9560K P=34;
i_vibr=1;
rel=1;
switch var
    case 1
    TT=readmatrix('T_50.csv');
    time_T_exp=TT(:,1);
    T_exp=TT(:,2);
    Tv_exp=TT(:,5);
    nm_exp=readmatrix('nm_50_1.csv');
    time_n_exp=nm_exp(:,1);
    nm_exp=nm_exp(:,2);
    time_p_exp=readmatrix('p.csv');
    p_exp=time_p_exp(:,2);
    time_p_exp=time_p_exp(:,1);
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
    case 7
    Test=readmatrix("6230K_Simulations.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,7);
    Tv_exp=Test(:,5);
    nm_exp=Test(:,9)*1e9;
    time_n_exp=time_T_exp;
    case 8
    Test=readmatrix("7940K_Simulations.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,7);
    Tv_exp=Test(:,5);
    nm_exp=Test(:,9)*1e9;
    time_n_exp=time_T_exp;
    case 9
    Test=readmatrix("9560K_Simulations.csv");
    time_T_exp=Test(:,1);
    T_exp=Test(:,7);
    Tv_exp=Test(:,5);
    nm_exp=Test(:,9)*1e9;
    time_n_exp=time_T_exp;
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

%nexttile 
%hold on
% semilogy(dat1(i_vibr,2,var,rel).time, ...
%     dat1(i_vibr,2,var,rel).nm_n./(dat1(i_vibr,2,var,rel).nm_n+dat1(i_vibr,2,var,rel).na+dat1(i_vibr,2,var,rel).nAr)...
%     ,'r-', 'LineWidth', 1.5, 'DisplayName', "n_m/n - U=D/6k ");
%  semilogy(dat1(i_vibr,3,var,rel).time, ...
%      dat1(i_vibr,3,var,rel).nm_n./(dat1(i_vibr,3,var,rel).nm_n+dat1(i_vibr,3,var,rel).na+dat1(i_vibr,3,var,rel).nAr)...
%      ,'b-', 'LineWidth', 1.5, 'DisplayName', "n_m/n - U=3T ");
%  semilogy(dat1(i_vibr,4,var,rel).time, ...
%      dat1(i_vibr,4,var,rel).nm_n./(dat1(i_vibr,4,var,rel).nm_n+dat1(i_vibr,4,var,rel).na+dat1(i_vibr,4,var,rel).nAr)...
%      ,'m-', 'LineWidth', 1.5, 'DisplayName', "n_m/n - U=inf ");
% legend('Location','ne');
% xlabel("t, \mus");
% ylabel("n_{O_2}/n_{sum}");
%hold off;
% grid minor;

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
%disp(Pressure);
%% 
% rel=2;
% i_vibr=1;
% for i=[1:9] 
%     tr(i,1)=dat(i_vibr,2,i).time(end);
% end
%%
% for i=[1:6]
%     P(i,1)=dat1(i_vibr, 2, i, 1).p(1);
%      P(i,2)=dat1(i_vibr, 3, i, 1).p(1);
%       P(i,3)=dat1(i_vibr, 4, i, 1).p(1);
%     P(i,4)=dat1(i_vibr, 2, i, 2).p(1);
%      P(i,5)=dat1(i_vibr, 3, i, 2).p(1);
%       P(i,6)=dat1(i_vibr, 4, i, 2).p(1);
% end

%%
%figure
%plot(dat(i_vibr, 4, 6).time, dat(i_vibr, 4, 6).nm./(dat(i_vibr, 4, 6).nm+dat(i_vibr, 4, 6).na+dat(i_vibr, 4, 6).nAr), 'k-.', 'LineWidth',1.5, "DisplayName","\chi_{O_2}, U=inf");
%plot(dat(i_vibr, 2, 6).time,log( dat(i_vibr, 2, 6).na./(dat(i_vibr, 2, 6).nm+dat(i_vibr, 2, 6).na+dat(i_vibr, 2, 6).nAr)), 'b--', 'LineWidth',1.5, "DisplayName","\chi_{O}, U=D/6k");
%plot(dat(i_vibr, 3, 6).time,log( dat(i_vibr,3, 6).na./(dat(i_vibr, 3, 6).nm+dat(i_vibr, 3, 6).na+dat(i_vibr, 3, 6).nAr)), 'm--', 'LineWidth',1.5, "DisplayName","\chi_{O}, U=3T");
%plot(dat(i_vibr, 4, 6).time,dat(i_vibr, 4, 6).na./(dat(i_vibr, 4, 6).nm+dat(i_vibr, 4, 6).na+dat(i_vibr, 4, 6).nAr), 'b--', 'LineWidth',1.5, "DisplayName","\chi_{O}, U=inf");
%plot(dat(i_vibr, 4, 6).time,dat(i_vibr, 4, 6).nAr./(dat(i_vibr, 4, 6).nm+dat(i_vibr, 4, 6).na+dat(i_vibr, 4, 6).nAr), 'r--', 'LineWidth',1.5, "DisplayName","\chi_{Ar}, U=inf");
% semilogy(dat(i_vibr, 4, 3).time, dat(i_vibr, 4, 3).nm./(dat(i_vibr, 4, 3).nm+dat(i_vibr, 4, 3).na+dat(i_vibr, 4, 3).nAr), 'k-.', dat(i_vibr, 4, 3).time,dat(i_vibr, 4, 3).na./(dat(i_vibr, 4, 3).nm+dat(i_vibr, 4, 3).na+dat(i_vibr, 4, 3).nAr), 'b--', dat(i_vibr, 4, 3).time,dat(i_vibr, 4, 3).nAr./(dat(i_vibr, 4, 3).nm+dat(i_vibr, 4, 3).na+dat(i_vibr, 4, 3).nAr), 'r--', 'LineWidth',1.5);
% xlabel("t, \mus");
% ylabel("\chi");
% legend("\chi_{O_2}, U=inf","\chi_{O}, U=inf","\chi_{Ar}, U=inf", 'Location','se');
% grid minor;

%%
% var=6;
% rel=1;
% i_vibr=1;
% figure
% hold on
% %title("Case " + info(var));
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "Tv - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "Tv - U=inf " );
% xlim([0 450]);
% legend('Location','se');
% xlabel("t, \mu s");
% ylabel("T_v, K");
% hold off
% grid minor

% figure
% semilogy( dat1(i_vibr,2,var,rel).time, ...
%   dat1(i_vibr,2,var,rel).nm_n./(dat1(i_vibr,2,var,rel).nm_n+dat1(i_vibr,2,var,rel).na+dat1(i_vibr,2,var,rel).nAr), 'b-', dat1(i_vibr,2,var,rel).time, ...
%   dat1(i_vibr,2,var,rel).na./(dat1(i_vibr,2,var,rel).nm_n+dat1(i_vibr,2,var,rel).na+dat1(i_vibr,2,var,rel).nAr), 'r-', dat1(i_vibr,3,var,rel).time, ...
%   dat1(i_vibr,3,var,rel).nm_n./(dat1(i_vibr,3,var,rel).nm_n+dat1(i_vibr,3,var,rel).na+dat1(i_vibr,3,var,rel).nAr), 'b--', dat1(i_vibr,3,var,rel).time, ...
%   dat1(i_vibr,3,var,rel).na./(dat1(i_vibr,3,var,rel).nm_n+dat1(i_vibr,3,var,rel).na+dat1(i_vibr,3,var,rel).nAr),'r--', dat1(i_vibr,4,var,rel).time, ...
%   dat1(i_vibr,4,var,rel).nm_n./(dat1(i_vibr,4,var,rel).nm_n+dat1(i_vibr,4,var,rel).na+dat1(i_vibr,4,var,rel).nAr), 'b-.', dat1(i_vibr,4,var,rel).time, ...
%   dat1(i_vibr,4,var,rel).na./(dat1(i_vibr,4,var,rel).nm_n+dat1(i_vibr,4,var,rel).na+dat1(i_vibr,4,var,rel).nAr), 'r-.','LineWidth', 1.5);
% 
% legend("\chi_{O2}, U=D/6k", "\chi_{O}, U=D/6k","\chi_{O2}, U=3T", "\chi_{O}, U=3T", "\chi_{O2}, U=\infty", "\chi_{O}, U=\infty", 'Location','sw');
% xlabel("t, \mus");
% ylabel("\chi");

%%
% rel = 1;
% figure 
% hold on
% p1=plot(time_T_exp, Tv_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', "Tv - exp case (" + num2str(var)+")");
% errorbar([time_T_exp(29) time_T_exp(43) time_T_exp(52) time_T_exp(63) time_T_exp(71)], [Tv_exp(29) Tv_exp(43) Tv_exp(52) Tv_exp(63) Tv_exp(71)],[Tv_exp(29) Tv_exp(43) Tv_exp(52) Tv_exp(63) Tv_exp(71)]*0.065, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',2, 'LineWidth', 1.5);
% %errorbar([time_T_exp(27) time_T_exp(39) time_T_exp(49) time_T_exp(56)], [Tv_exp(27) Tv_exp(39) Tv_exp(49) Tv_exp(56)],[Tv_exp(27) Tv_exp(39) Tv_exp(49) Tv_exp(56)]*0.038, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',2, 'LineWidth', 1.5);
% p2=plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', 1.5, 'DisplayName', "U=D/6k SSH" );
% p3=plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).Tv,'b--', 'LineWidth', 1.5, 'DisplayName', "U=3T SSH" );
% p4=plot(dat1(2,3,var,rel).time, dat1(2,3,var,rel).Tv,'b-', 'LineWidth', 1.5, 'DisplayName', "U=3T FHO" )
% p5=plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).Tv, 'm-','LineWidth', 1.5, 'DisplayName', "U=\infty SSH" );
% xlim([0 50]);
% legend([p1, p2, p3, p4, p5],'Location','se');
% xlabel("t, \mu s");
% ylabel("T_v, K");
% hold off
% grid minor

%%
% rel=2;
% figure 
% hold on
% p1=plot(time_n_exp, nm_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', "n_m - exp case " + num2str(var));
% errorbar([time_n_exp(4) time_n_exp(6) time_n_exp(7)], [nm_exp(4) nm_exp(6) nm_exp(7)],[nm_exp(4)*0.024 nm_exp(6)*0.027 nm_exp(7)*0.036], 'sqk', 'MarkerFaceColor', 'k','MarkerSize',2, 'LineWidth', 1.5);
% %errorbar([time_n_exp(4) time_n_exp(5) time_n_exp(6)], [nm_exp(4) nm_exp(5) nm_exp(6)],[nm_exp(4)*0.029 nm_exp(6)*0.042 nm_exp(7)*0.06], 'sqk', 'MarkerFaceColor', 'k','MarkerSize',2, 'LineWidth', 1.5);
% p2=plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).nm_n*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "U=D/6k SSH" );
% p3=plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).nm_n*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "U=3T SSH" );
% p4=plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).nm_n*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "U=inf SSH" );
% legend([p1 p2 p3 p4],'Location','ne');
% xlabel("t, \mu s");
% ylabel("n_{O_2}, mmol/m^3");
% xlim([0 50]);
% hold off
% grid minor
