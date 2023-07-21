%%Comparison with experimental data and plots
%Streicher's 2022 experiment with NO\O\N\O2\N2\Ar mixtures
clearvars;
load('..\data\NO_N2 Streicher experiment\output_ReflSW_NO_N2.mat');
load('..\data\NO_N2 Streicher experiment\NO_N2_Streicher22_experiment.mat');
info=["2% NO/N2 No.1", "2% NO/N2 No.2", "2% NO/N2 No.3" ,"0.4% NO/N2 No.4" ,"0.4% NO/N2 No.5", "0.4% NO/N2 No.6",...
    "2% NO/N2/Ar No.7","2% NO/N2/Ar No.8","2% NO/N2/Ar No.9","0.4% NO/N2/Ar No.10","0.4% NO/N2/Ar No.11", "0.4% NO/N2/Ar No.12"];

Na=6.02214076e23;
%Calculated data in case with exchange reactions and without vibrational
%activation of the reaction product
%dat1=data_behindRSW_withexch_VDOP0; 
%dat=data_betweenSWs_withexch_VDOP0;

%Calculated data in case with exchange reactions and with vibrational
%activation of the reaction product
dat1=data_behindRSW_withexch_VDOP1; 
dat=data_betweenSWs_withexch_VDOP1;

%Calculated data in case without exchange reactions
%dat1=data_behindRSW_withoutexch;
%dat=data_betweenSWs_withoutexch;

for var=5
%testcases     
%var: %1 - 2-02 T=3560 P=0.561;  2 - 2-14 T=5460 P=0.325; 3 - 2-32 T=7070
%P=0.119; 4 - 2-38 T=8730 P=0.137
% 5 - 1-01 T=2360 P=0.757; 6 - 1-04 T=3470 P=0.584; 7 - 1-08 T=5580K P=0.252;
% 8 - 1-12 T=7090K P=0.274; 9 - 04-02 T=2180K P=0.952; 10 - 04-11 T=3470K 
%P=0.668;  11 - 04-22 T=5650K P=0.545;
i_vibr=2; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on
time_T_err=[];
T_err=[];
err_T=[];
time_n_err=[];
n_err=[];
err_n=[];
time_n0_err=[];
n0_err=[];
err_n0=[];
%initialization of experimental data for temperature
j=1;
while (data_experiment(var).T(j,3)~=0)
j=j+1;
end
time_T_exp=data_experiment(var).T(1:j-1,1);
T_exp=data_experiment(var).T(1:j-1,2);
[time_T_exp, I]=sort(time_T_exp);
T_exp=T_exp(I);
tlim_T=max(time_T_exp);
for i=j+1:length(data_experiment(var).T(:,3))-1
    if (data_experiment(var).T(i,3)~=-1)&&(data_experiment(var).T(i-1,3)==data_experiment(var).T(i,3))&&(data_experiment(var).T(i+1,3)==data_experiment(var).T(i,3))
        time_T_err=[time_T_err, (data_experiment(var).T(i,1)+data_experiment(var).T(i-1,1)+data_experiment(var).T(i+1,1))/3];
        T_err=[T_err, data_experiment(var).T(i-1,2)];
        err_T=[err_T, (data_experiment(var).T(i,2) - data_experiment(var).T(i+1,2))/2];
    end
end

time_TvNO_exp=data_experiment(var).TvNO(:,1);
TvNO_exp=data_experiment(var).TvNO(:,2);
[time_TvNO_exp, I]=sort(time_TvNO_exp);
TvNO_exp=TvNO_exp(I);

time_TvN2_exp=data_experiment(var).TvN2(:,1);
TvN2_exp=data_experiment(var).TvN2(:,2);
[time_TvN2_exp, I]=sort(time_TvN2_exp);
TvN2_exp=TvN2_exp(I);

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

j=1;
while (data_experiment(var).n0(j,3)~=0)
j=j+1;
end
time_n0_exp=data_experiment(var).n0(1:j-1,1);
n0_exp=data_experiment(var).n0(1:j-1,2);
[time_n0_exp, I]=sort(time_n0_exp);
n0_exp=n0_exp(I);
for i=j+1:length(data_experiment(var).n0(:,3))-1
    if (data_experiment(var).n0(i,3)~=-1)&&(data_experiment(var).n0(i-1,3)==data_experiment(var).n0(i,3))&&(data_experiment(var).n0(i+1,3)==data_experiment(var).n0(i,3))
        time_n0_err=[time_n0_err, (data_experiment(var).n0(i,1)+data_experiment(var).n0(i-1,1)+data_experiment(var).n0(i+1,1))/3];
        n0_err=[n0_err, data_experiment(var).n0(i-1,2)];
        err_n0=[err_n0, (data_experiment(var).n0(i,2) - data_experiment(var).n0(i+1,2))/2];
    end
end

%err=0.045;
   

% figure("Position", [0, 0, 1400, 800])
% t=tiledlayout(2, 3, "TileSpacing", "compact");
% title(t, "Case " + info(var));
% nexttile
% hold on
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
% plot(time_T_exp, T_exp, 'k-', 'LineWidth', 2, 'DisplayName', "T - experiment ");
% errorbar(time_T_err, T_err, err_T, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1)
% legend('Location','best');
% xlim([-10 tlim_T]);
% ylim([0 max(T_exp)+500]);
% xlabel("t, \mu s");
% ylabel("T, K");
% hold off
% grid minor
% nexttile
% hold on
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).TvNO,'r-', 'LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).TvNO,'b-', 'LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).TvNO, 'm-','LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=inf " );
% plot(time_TvNO_exp, TvNO_exp, 'k-', 'LineWidth', 2, 'DisplayName', "Tv - experiment");
% xlim([-10 tlim_T]);
% %ylim([0 6000]);
% legend('Location','best');
% xlabel("t, \mu s");
% ylabel("T_v^{NO}, K");
% hold off
% grid minor
% nexttile
% hold on
% plot(time_n_exp, n_exp, 'k-', 'LineWidth', 2, 'DisplayName', "n_{NO} - experiment" + num2str(var));
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).nNO*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n_{NO} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).nNO*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n_{NO} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).nNO*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n_{NO} - U=inf " );
% errorbar(time_n_err, n_err, err_n, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
% xlim([-10 tlim_n]);
% ylim([0 max(dat1(i_vibr,2,var,rel).nNO*1e3)+4]);
% legend('Location','best');
% xlabel("t, \mu s");
% ylabel("n_{NO}, mmol/m^3");
% hold off
% grid minor
% nexttile
% hold on
% plot(time_n0_exp, n0_exp, 'k-', 'LineWidth', 2, 'DisplayName', "n^0_{NO} - experiment" + num2str(var));
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).ni_NO(:,1)/Na*1e3,'r-', 'LineWidth', 1.5, 'DisplayName', "n^0_{NO} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).ni_NO(:,1)/Na*1e3,'b-', 'LineWidth', 1.5, 'DisplayName', "n^0_{NO} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).ni_NO(:,1)/Na*1e3, 'm-','LineWidth', 1.5, 'DisplayName', "n^0_{NO} - U=inf " );
% errorbar(time_n0_err, n0_err, err_n0, 'sqk', 'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
% xlim([-10 tlim_n]);
% ylim([0 max(dat1(i_vibr,2,var,rel).nNO*1e3)+4]);
% legend('Location','best');
% xlabel("t, \mu s");
% ylabel("n^0_{NO}, mmol/m^3");
% hold off
% grid minor
% nexttile
% hold on
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).TvN2,'r-', 'LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).TvN2,'b-', 'LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).TvN2, 'm-','LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=inf " );
% plot(time_TvN2_exp, TvN2_exp, 'k-', 'LineWidth', 2, 'DisplayName', "Tv - experiment");
% xlim([-10 tlim_T]);
% %ylim([0 6000]);
% legend('Location','best');
% xlabel("t, \mu s");
% ylabel("T_v^{N_2}, K");
% hold off
% grid minor
% nexttile
% hold on
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).T,'r-', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).T,'b-', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).T, 'm-','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).TvNO,'r--', 'LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).TvNO,'b--', 'LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).TvNO, 'm--','LineWidth', 1.5, 'DisplayName', "T_v^{NO} - U=inf " );
% plot(dat1(i_vibr,2,var,rel).time, dat1(i_vibr,2,var,rel).TvN2,'r-.', 'LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=D/6k " );
% plot(dat1(i_vibr,3,var,rel).time, dat1(i_vibr,3,var,rel).TvN2,'b-.', 'LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=3T " );
% plot(dat1(i_vibr,4,var,rel).time, dat1(i_vibr,4,var,rel).TvN2, 'm-.','LineWidth', 1.5, 'DisplayName', "T_v^{N2} - U=inf " );
% xlim([-10 600]);
% %ylim([0 6000]);
% legend('Location','best');
% xlabel("t, \mu s");
% ylabel("T, K");
% hold off
% grid minor

addpath('../src/');
addpath('../data/')
load('particles.mat', "NO", "N", "O", "Ar", "O2", "N2");
k=1.380649e-23; Na=6.02214076e23;
NO.num_elex_levels=1;
N2.num_elex_levels=1;
O2.num_elex_levels=1;
i_U=2;
[~,j1]=min(abs(dat1(i_vibr, i_U, var, rel).time - tlim_n/3));
[~,j2]=min(abs(dat1(i_vibr, i_U, var, rel).time - 2*tlim_n/3));
[~,j3]=min(abs(dat1(i_vibr, i_U, var, rel).time - tlim_n));

figure("Position", [200, 100, 900, 800])
t=tiledlayout(2, 2, "TileSpacing", "compact");
title(t, "Case " + info(var) + ". Max T_v^{NO}(t)=" + num2str(max((dat1(i_vibr, i_U, var, rel).TvNO))));
nexttile
n_bolz=density_f_exc(dat1(i_vibr, i_U, var, rel).TvNO(3), dat1(i_vibr, i_U, var, rel).nNO(3)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz, 'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withoutexch(i_vibr, i_U, var, rel).ni_NO(3, :) ,'r-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP1(i_vibr, i_U, var, rel).ni_NO(3, :) , 'm-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP0(i_vibr, i_U, var, rel).ni_NO(3, :) , 'g-', 'LineWidth',2.0);
xlabel('NO vibr. level, i');
ylabel('n^i_{NO}, м^{-3}');
lgd=legend('Распределение по Больцману','Распределение без обм. р-й','Распределение с обм. р-ми VDOP=1','Распределение с обм. р-ми VDOP=0','Location','south');
xlim([0 NO.num_vibr_levels(1)-1]);
title("t=0 \mus" + "  T_v^{NO}(t)="+num2str(dat1(i_vibr, i_U, var, rel).TvNO(3)));
hold off
grid minor

nexttile
n_bolz=density_f_exc(dat1(i_vibr, i_U, var, rel).TvNO(j1), dat1(i_vibr, i_U, var, rel).nNO(j1)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz, 'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withoutexch(i_vibr, i_U, var, rel).ni_NO(j1, :) ,'r-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP1(i_vibr, i_U, var, rel).ni_NO(j1, :) , 'm-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP0(i_vibr, i_U, var, rel).ni_NO(j1, :) , 'g-', 'LineWidth',2.0);
xlabel('NO vibr. level, i');
ylabel('n^i_{NO}, м^{-3}');
%legend('Распределение по Больцману','Распределение без обм. р-й','Распределение с обм. р-ми VDOP=1','Распределение с обм. р-ми VDOP=0','Location','south');
xlim([0 NO.num_vibr_levels(1)-1]);
title("t="+num2str(round(tlim_n/3))+" \mus"  + "  T_v^{NO}(t)="+num2str(dat1(i_vibr, i_U, var, rel).TvNO(j1)))
hold off
grid minor

nexttile
n_bolz=density_f_exc(dat1(i_vibr, i_U, var, rel).TvNO(j2), dat1(i_vibr, i_U, var, rel).nNO(j2)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz, 'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withoutexch(i_vibr, i_U, var, rel).ni_NO(j2, :) ,'r-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP1(i_vibr, i_U, var, rel).ni_NO(j2, :) , 'm-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP0(i_vibr, i_U, var, rel).ni_NO(j2, :) , 'g-', 'LineWidth',2.0);
xlabel('NO vibr. level, i');
ylabel('n^i_{NO}, м^{-3}');
%legend('Распределение по Больцману','Распределение без обм. р-й','Распределение с обм. р-ми VDOP=1','Распределение с обм. р-ми VDOP=0','Location','south');
xlim([0 NO.num_vibr_levels(1)-1]);
title("t="+num2str(round(2*tlim_n/3))+" \mus" + "  T_v^{NO}(t)="+num2str(dat1(i_vibr, i_U, var, rel).TvNO(j2)))
hold off
grid minor

nexttile
n_bolz=density_f_exc(dat1(i_vibr, i_U, var, rel).TvNO(j3), dat1(i_vibr, i_U, var, rel).nNO(j3)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz, 'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withoutexch(i_vibr, i_U, var, rel).ni_NO(j3, :) ,'r-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP1(i_vibr, i_U, var, rel).ni_NO(j3, :) , 'm-','LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, data_behindRSW_withexch_VDOP0(i_vibr, i_U, var, rel).ni_NO(j3, :) , 'g-', 'LineWidth',2.0);
xlabel('NO vibr. level, i');
ylabel('n^i_{NO}, м^{-3}');
%legend('Распределение по Больцману','Распределение без обм. р-й','Распределение с обм. р-ми VDOP=1','Распределение с обм. р-ми VDOP=0','Location','south');
xlim([0 NO.num_vibr_levels(1)-1]);
title("t="+num2str(round(tlim_n))+" \mus" + "  T_v^{NO}(t)="+num2str(dat1(i_vibr, i_U, var, rel).TvNO(j3)))
hold off
grid minor

lgd.Layout.Tile = 'south';
rmpath('../src/')
rmpath('../data/')
end