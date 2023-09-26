%%Comparison with experimental data and plots
%Streicher's 2022 experiment with NO\O\N\O2\N2\Ar mixtures


clearvars;
%Here you need to load an array with calculated data
load('..\data\NO Streicher experiment\NO_behind_ReflSW.mat');
load('..\data\NO Streicher experiment\NO_Streicher22_experiment.mat');
info=["No.1: NO - 2%; Ar - 98%", "2% No.2", "2% No.3" ,"2% No.4" ,...
    "1% No.5", "1% No.6", "No.7: NO - 1%; Ar - 99%","1% No.8",...
    "0.4% No.9","0.4% No.10","0.4% No.11"];
Na=6.02214076e23;


%Initialization of calculated data. data_main - the main model, 
% in our case it was a model with included exchange reactions 
% according to the Kunova model. data_minor - is data for comparison 
% with the data_main, in our case these were models without taking into 
% account exchange reactions or with averaging of the Kunova model
data_main=dat1(:,:,:,:,2);
data_minor=dat1(:,:,:,:,1);

%indicators for plots. True if you want to plot a graph,
%False if you dont want to.
%MACRO_plot - plot with macro parameters: T, TvNO, n_NO, n_NO_0
%NNO_plot - plot with comparison of different models for n_NO
%VDF_NO_plot - plot VDF for NO
%TVerror_plot - plot deviation for Tv_NO
%Nerror_plot - plot deviation for n_NO

MACRO_plot=true;
NNO_plot=true;
VDF_NO_plot=true;
TVerror_plot=true;
Nerror_plot=true;

for var=9
%testcases     
%var: %1 - 2-02 T=3560 P=0.561;  2 - 2-14 T=5460 P=0.325; 3 - 2-32 T=7070
%P=0.119; 4 - 2-38 T=8730 P=0.137
% 5 - 1-01 T=2360 P=0.757; 6 - 1-04 T=3470 P=0.584; 7 - 1-08 T=5580K P=0.252;
% 8 - 1-12 T=7090K P=0.274; 9 - 04-02 T=2180K P=0.952; 10 - 04-11 T=3470K 
%P=0.668;  11 - 04-22 T=5650K P=0.545;
i_vibr=2; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on

%initialization of experimental data for temperature
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
    if (data_experiment(var).T(i,3)~=-1)&&...
        (data_experiment(var).T(i-1,3)==data_experiment(var).T(i,3))&&...
        (data_experiment(var).T(i+1,3)==data_experiment(var).T(i,3))
        time_T_err=[time_T_err, (data_experiment(var).T(i,1)+ ...
           data_experiment(var).T(i-1,1)+data_experiment(var).T(i+1,1))/3];
        T_err=[T_err, data_experiment(var).T(i-1,2)];
        err_T=[err_T, (data_experiment(var).T(i,2) - ...
            data_experiment(var).T(i+1,2))/2];
    end
end

time_TvNO_exp=data_experiment(var).TvNO(:,1);
TvNO_exp=data_experiment(var).TvNO(:,2);
[time_TvNO_exp, I]=sort(time_TvNO_exp);
TvNO_exp=TvNO_exp(I);

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
    if (data_experiment(var).n(i,3)~=-1)&&...
        (data_experiment(var).n(i-1,3)==data_experiment(var).n(i,3))&&...
        (data_experiment(var).n(i+1,3)==data_experiment(var).n(i,3))
        time_n_err=[time_n_err, (data_experiment(var).n(i,1)+ ...
           data_experiment(var).n(i-1,1)+data_experiment(var).n(i+1,1))/3];
        n_err=[n_err, data_experiment(var).n(i-1,2)];
        err_n=[err_n, (data_experiment(var).n(i,2) - ...
            data_experiment(var).n(i+1,2))/2];
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
    if (data_experiment(var).n0(i,3)~=-1)&&...
        (data_experiment(var).n0(i-1,3)==data_experiment(var).n0(i,3))&&...
        (data_experiment(var).n0(i+1,3)==data_experiment(var).n0(i,3))
        time_n0_err=[time_n0_err, (data_experiment(var).n0(i,1)+ ...
         data_experiment(var).n0(i-1,1)+data_experiment(var).n0(i+1,1))/3];
        n0_err=[n0_err, data_experiment(var).n0(i-1,2)];
        err_n0=[err_n0, (data_experiment(var).n0(i,2) - ...
            data_experiment(var).n0(i+1,2))/2];
    end
end
%End of experimental data handling

if MACRO_plot
%Macro parameters plot
figure("Position", [0, 0, 900, 800])
t=tiledlayout(2, 2, "TileSpacing", "compact");
%title(t, "Case " + info(var) + "; vibr. model: FHO",...
% 'FontName', 'Times New Roman');

nexttile
hold on
set(gca, 'FontName', 'Times New Roman');
p1=plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).T,'r--', 'LineWidth', ...
    1.5, 'DisplayName', "\it T - U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).T,'b--', 'LineWidth', ...
    1.5, 'DisplayName', "\it T - U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).T, 'm--','LineWidth', ...
    1.5, 'DisplayName', "\it T - U=\infty " );
errorbar(time_T_err, T_err, err_T, 'sqk', ...
    'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
time_T_err=[0 time_T_err tlim_T];
ppp=spline(time_T_exp, T_exp, time_T_err);
p4=plot(time_T_err, ppp, 'k-', 'LineWidth', 2, ...
    'DisplayName', "\it T - \rm experiment ");
legend([p1 p2 p3 p4],'Location','best');
xlim([-10 tlim_T]);
ylim([min(T_exp)-300 max(T_exp)+300]);
xlabel("\it t\rm, \mus");
ylabel("\it T\rm, K");
hold off
grid minor

nexttile
hold on
set(gca, 'FontName', 'Times New Roman');
plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).TvNO,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it T}_v^{NO} - \it U=D/6k " );
plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).TvNO,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it T}_v^{NO} - \it U=3T " );
plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).TvNO, 'm-','LineWidth', ...
    1.5, 'DisplayName', "{\it T}_v^{NO} - \it U=\infty " );
plot(time_TvNO_exp, TvNO_exp, 'k-', 'LineWidth', ...
    2, 'DisplayName', "{\it T}_v^{NO} - experiment");
xlim([-10 tlim_T]);
%ylim([0 6000]);
legend('Location','best');
xlabel("\it t\rm, \mus");
ylabel("{\it T}_v^{NO}, K");
hold off
grid minor

nexttile
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(time_n_exp, n_exp, 'k-', 'LineWidth', 2, ...
    'DisplayName', "\it n_{\rm NO} \rm - experiment" );
p2=plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nNO*1e3,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm NO} - U=D/6k " );
p3=plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nNO*1e3,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm NO} - U=3T " );
p4=plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nNO*1e3, 'm-','LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm NO} - U=\infty " );
errorbar(time_n_err, n_err, err_n, 'sqk', 'MarkerFaceColor', ...
    'k','MarkerSize',1, 'LineWidth', 1);
xlim([-10 tlim_n]);
ylim([0 max(data_main(i_vibr,2,var,rel).nNO*1e3)+4]);
legend([p1 p2 p3 p4],'Location','best');
xlabel("\it t\rm, \mus");
ylabel("{\it n}_{NO}, mmol/m^3");
hold off
grid minor

nexttile
hold on
set(gca, 'FontName', 'Times New Roman');
p1=plot(time_n0_exp, n0_exp, 'k-', 'LineWidth', 2, 'DisplayName', ...
    "{\it n}^0_{NO} - experiment");
p2=plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).ni_NO(:,1)/Na*1e3,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it n}^0_{NO} - \it U=D/6k " );
p3=plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).ni_NO(:,1)/Na*1e3,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it n}^0_{NO} - \it U=3T " );
p4=plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).ni_NO(:,1)/Na*1e3, 'm-','LineWidth', ...
    1.5, 'DisplayName', "{\it n}^0_{NO} - \it U=\infty " );
errorbar(time_n0_err, n0_err, err_n0, 'sqk', 'MarkerFaceColor', ...
    'k','MarkerSize',1, 'LineWidth', 1);
xlim([-10 tlim_n]);
legend([p1 p2 p3 p4],'Location','best');
xlabel("\it t\rm, \mus");
ylabel("{\it n}^0_{NO}, mmol/m^3");
hold off
grid minor
end

if NNO_plot
figure
hold on
p1=plot(time_n_exp, n_exp, 'k-', 'LineWidth', 2, ...
    'DisplayName', "{\it n}_{NO} - experiment");
g1=spline(data_main(i_vibr,2,var,rel).time, ...
   data_main(i_vibr,2,var,rel).nNO*1e3, data_minor(i_vibr,2,var,rel).time);
g2=spline(data_main(i_vibr,3,var,rel).time, ...
   data_main(i_vibr,3,var,rel).nNO*1e3, data_minor(i_vibr,3,var,rel).time);
g3=spline(data_main(i_vibr,4,var,rel).time, ...
   data_main(i_vibr,4,var,rel).nNO*1e3, data_minor(i_vibr,4,var,rel).time);
p2=plot(data_minor(i_vibr,2,var,rel).time, g1,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it n}_{NO} -\it U=D/6k VDOP1" );
p3=plot(data_minor(i_vibr,3,var,rel).time, g2,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "{\it n}_{NO} -\it U=3T VDOP1" );
p4=plot(data_minor(i_vibr,4,var,rel).time, g3, 'm-','LineWidth', ...
    1.5, 'DisplayName', "{\it n}_{NO} -\it U=\infty VDOP1" );
dev1=max(abs(g1-data_minor(i_vibr,2,var,rel).nNO*1e3)./g1)
dev2=max(abs(g2-data_minor(i_vibr,3,var,rel).nNO*1e3)./g2)
dev3=max(abs(g3-data_minor(i_vibr,4,var,rel).nNO*1e3)./g3)
dev4=mean(abs(g1-data_minor(i_vibr,2,var,rel).nNO*1e3)./g1)
dev5=mean(abs(g2-data_minor(i_vibr,3,var,rel).nNO*1e3)./g2)
dev6=mean(abs(g3-data_minor(i_vibr,4,var,rel).nNO*1e3)./g3)
p5=plot(data_minor(i_vibr,2,var,rel).time, ...
    data_minor(i_vibr,2,var,rel).nNO*1e3,'r--', 'LineWidth', 1.5, ...
    'DisplayName', "{\it n}_{NO} -\it U=D/6k VDOP0" );
p6=plot(data_minor(i_vibr,3,var,rel).time, ...
    data_minor(i_vibr,3,var,rel).nNO*1e3,'b--', 'LineWidth', ...
    1.5, 'DisplayName', "{\it n}_{NO} -\it U=3T VDOP0" );
p7=plot(data_minor(i_vibr,4,var,rel).time, ...
    data_minor(i_vibr,4,var,rel).nNO*1e3, 'm--','LineWidth', ...
    1.5, 'DisplayName', "{\it n}_{NO} -\it U=\infty VDOP0" );
errorbar(time_n_err, n_err, err_n, 'sqk', 'MarkerFaceColor', ...
    'k','MarkerSize',1, 'LineWidth', 1);
xlim([-10 tlim_n]);
ylim([0 max(data_main(i_vibr,2,var,rel).nNO*1e3)+4]);
legend([p1 p2 p3 p4],'Location','best');
xlabel("{\it t}, \mus");
ylabel("{\it n}_{NO}, mmol/m^3");
hold off
grid minor
end

addpath('../src/');
addpath('../data/')
load('particles.mat', "NO", "N", "O", "Ar", "O2", "N2");
NO.num_elex_levels=1;
N2.num_elex_levels=1;
O2.num_elex_levels=1;
k=1.380649e-23; Na=6.02214076e23;
i_U=2; %Choose the dissotiation model

%total number density of the mixture in case without exch. reactions
nSUMwoexch=(data_minor(i_vibr, i_U, var, rel).nNO + ...
    data_minor(i_vibr, i_U, var, rel).nN2 + ...
    data_minor(i_vibr, i_U, var, rel).nO2 + ...
    data_minor(i_vibr, i_U, var, rel).nN + ...
    data_minor(i_vibr, i_U, var, rel).nO + ...
    data_minor(i_vibr, i_U, var, rel).nAr)*Na;

%total number density of the mixture in case with exch. reactions
nSUMwexch=(data_main(i_vibr, i_U, var, rel).nNO + ...
    data_main(i_vibr, i_U, var, rel).nN2 + ...
    data_main(i_vibr, i_U, var, rel).nO2 + ...
    data_main(i_vibr, i_U, var, rel).nN + ...
    data_main(i_vibr, i_U, var, rel).nO + ...
    data_main(i_vibr, i_U, var, rel).nAr)*Na;

[~,j1]=min(abs(data_main(i_vibr, i_U, var, rel).time - tlim_n/3));
[~,j2]=min(abs(data_main(i_vibr, i_U, var, rel).time - 2*tlim_n/3));
[~,j3]=min(abs(data_main(i_vibr, i_U, var, rel).time - tlim_n));

[~,j11]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - tlim_n/3));
[~,j22]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - 2*tlim_n/3));
[~,j33]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - tlim_n));
if VDF_NO_plot
ylim_low=1e-15;
figure("Position", [200, 100, 900, 800])
t=tiledlayout(2, 2, "TileSpacing", "compact");
% title(t, "Case " + info(var) + ". Max T_v^{NO}(t)=" + ...
%     num2str(max((data_main(i_vibr, i_U, var, rel).TvNO))));
nexttile %t=0 
n_bolz=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(1),...
    data_main(i_vibr, i_U, var, rel).nNO(1)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz./nSUMwexch(1), 'k-',...
    'LineWidth',2.0, "DisplayName","Boltzmann");
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:NO.num_vibr_levels(1)-1,...
    data_minor(i_vibr, i_U, var, rel).ni_NO(1, :)./...
    nSUMwoexch(1) ,'-','color',[0 0.6 0],'LineWidth',2.0,...
    "DisplayName","Without exch. react.");
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_main(i_vibr, i_U, var, rel).ni_NO(1, :)./...
    nSUMwexch(1) , 'm-','LineWidth',2.0, "DisplayName",...
    "With exch. react. and vibr. activation of react. prod. (NO)");
xlabel('NO vibr. level,\it i');
ylabel('\it n_{\rm NO\iti}/\itn');
lgd=legend('Location','south');
xlim([0 NO.num_vibr_levels(1)-1]);
ylim([ylim_low*1e-15 1])
text(3,ylim_low*1e-10,"{\it t} = 0 \mus," + " {\it T}_v^{NO}({\itt}) = "...
    +num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(1)))+" K",...
    'FontSize',12, 'FontName', 'Times New Roman');
hold off
box off
grid minor

nexttile %t=tlim/3
n_bolz=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j1),...
    data_main(i_vibr, i_U, var, rel).nNO(j1)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz./nSUMwexch(j1), 'k-', ...
    'LineWidth',2.0);
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:NO.num_vibr_levels(1)-1,...
    data_minor(i_vibr, 2, var, rel).ni_NO(j11, :)./...
    nSUMwoexch(j11) ,'-','color',[0 0.6 0],'LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1,...
    data_main(i_vibr, 2, var, rel).ni_NO(j1, :)./...
    nSUMwexch(j1) , 'm-','LineWidth',2.0);
xlabel('NO vibr. level,\it i');
ylabel('\it n_{\rm NO\iti}/\itn');
xlim([0 NO.num_vibr_levels(1)-1]);
ylim([ylim_low 1])
text(3,ylim_low*1e2,"{\it t} = "+num2str(round(tlim_n/3))+" \mus,"  + ...
    " {\it T}_v^{NO}({\itt}) = "+ ...
    num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(j1)))+" K", ...
    'FontSize',12, 'FontName', 'Times New Roman')
hold off
box off
grid minor

nexttile %t=2*tlim/3
n_bolz=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j2), ...
    data_main(i_vibr, i_U, var, rel).nNO(j2)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz./nSUMwexch(j2), ...
    'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_minor(i_vibr, i_U, var, rel).ni_NO(j22, :)./ ...
    nSUMwoexch(j22) ,'-','color',[0 0.6 0],'LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_main(i_vibr, i_U, var, rel).ni_NO(j2, :)./ ...
    nSUMwoexch(j2) , 'm-','LineWidth',2.0);
xlabel('NO vibr. level,\it i');
ylabel('\it n_{\rm NO\iti}/\itn');
xlim([0 NO.num_vibr_levels(1)-1]);
ylim([ylim_low 1])
text(3,ylim_low*1e2,"{\it t} = "+num2str(round(2*tlim_n/3))+" \mus," ...
    + " {\it T}_v^{NO}({\itt}) = "+ ...
    num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(j2)))+" K", ...
    'FontSize',12, 'FontName', 'Times New Roman')
hold off
box off
grid minor

nexttile %t=tlim
n_bolz=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j3), ...
    data_main(i_vibr, i_U, var, rel).nNO(j3)*Na, NO);
semilogy(0:NO.num_vibr_levels(1)-1, n_bolz./nSUMwoexch(j3), 'k-', ...
    'LineWidth',2.0);
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_minor(i_vibr, i_U, var, rel).ni_NO(j33, :)./ ...
    nSUMwoexch(j33) ,'-','color',[0 0.6 0],'LineWidth',2.0);
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_main(i_vibr, i_U, var, rel).ni_NO(j3, :)./ ...
    nSUMwoexch(j3) , 'm-','LineWidth',2.0);
xlabel('NO vibr. level,\it i', 'FontName','Times New Roman');
ylabel('\it n_{\rm NO\iti}/\itn');
xlim([0 NO.num_vibr_levels(1)-1]);
ylim([ylim_low 1])
text(3, ylim_low*1e2,"{\it t} = "+num2str(round(tlim_n))+" \mus," ...
    + " {\it T}_v^{NO}({\itt}) = "+ ...
    num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(j3)))+" K", ...
    'FontSize',12, 'FontName', 'Times New Roman')
hold off
box off
grid minor
lgd.Layout.Tile = 'south';
end

rmpath('../src/')
rmpath('../data/')

end
%%

TEMP=[3560 5460 7070 8730 2360 3470 5580 7090 2180 3470 5650];
exp_err=[0.016 0.15 0.55 0.65 0.012 0.03 0.25 0.55 0.024 0.08 0.44];

if TVerror_plot
figure
hold on
for var=1:11

%separation of data and errors
time_TvNO_exp=data_experiment(var).TvNO(:,1);
TvNO_exp=data_experiment(var).TvNO(:,2);
[time_TvNO_exp, I]=sort(time_TvNO_exp);
TvNO_exp=TvNO_exp(I);
k=1;

%throwing out points with negative time
for i=1:length(TvNO_exp)
    if time_TvNO_exp(i)<0
        k=i;
    end
end
time_TvNO_exp=time_TvNO_exp(k+1:end);
TvNO_exp=TvNO_exp(k+1:end);
%throwing out points with the same abscissa
i=2;
while i~=length(time_TvNO_exp)
    if time_TvNO_exp(i)==time_TvNO_exp(i-1)
        time_TvNO_exp=[time_TvNO_exp(1:i-1),
            time_TvNO_exp(i+1:end)];
        TvNO_exp=[TvNO_exp(1:i-1),
            TvNO_exp(i+1:end)];
    end
    i=i+1;
end
%interpolation of calculated data for experimental time
i_vibr=1; rel=2;
s1=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH D/6k
s2=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH 3T
s3=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH inf
i_vibr=2;
s4=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO D/6k
s5=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO 3T
s6=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).TvNO, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO inf

%Uncomment if you want to account cases without relaxation
%dont forgot to fix a legend
% rel=1;
% i_vibr=1;
% s7=spline(data_main(i_vibr,2,var,rel).time, ...
%     data_main(i_vibr,2,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH D/6k
% s8=spline(data_main(i_vibr,3,var,rel).time, ...
%     data_main(i_vibr,3,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH 3T
% s9=spline(data_main(i_vibr,4,var,rel).time, ...
%     data_main(i_vibr,4,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %SSH inf
% i_vibr=2;
% s10=spline(data_main(i_vibr,2,var,rel).time, ...
%     data_main(i_vibr,2,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO D/6k
% s11=spline(data_main(i_vibr,3,var,rel).time, ...
%     data_main(i_vibr,3,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO 3T
% s12=spline(data_main(i_vibr,4,var,rel).time, ...
%     data_main(i_vibr,4,var,rel).TvNO, ...
%     min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %FHO inf


s_exp=makima(time_TvNO_exp, TvNO_exp, ...
    min(time_TvNO_exp):0.1:max(time_TvNO_exp)); %EXP

%calculating of mean deviation
err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
%Uncomment if you want to account cases without relaxation
%dont forgot to fix a legend
% err7=mean(abs(s7-s_exp)./s_exp);
% err8=mean(abs(s8-s_exp)./s_exp);
% err9=mean(abs(s9-s_exp)./s_exp);
% err10=mean(abs(s10-s_exp)./s_exp);
% err11=mean(abs(s11-s_exp)./s_exp);
% err12=mean(abs(s12-s_exp)./s_exp);

plot(TEMP(var), err1, 'ok');
plot(TEMP(var), err4, 'o','color', [0 0.6 0]);
%plot(TEMP(var), err7, 'ob');
%plot(TEMP(var), err10,'o','color', [0.9 0 0]);
plot(TEMP(var), err2, 'sk');
plot(TEMP(var), err5, 's','color',[0 0.6 0]);
%plot(TEMP(var), err8, 'sb');
%plot(TEMP(var), err11, 's','color', [0.9 0 0]);
plot(TEMP(var), err3, 'dk');
plot(TEMP(var), err6,'d','color',[0 0.6 0]);
%plot(TEMP(var), err9, 'db');
%plot(TEMP(var), err12, 'd','color', [0.9 0 0]);
end
hold off
legend("SSH, D/6k", "FHO, D/6k"...
 , "SSH, 3T", "FHO, 3T"...
 , "SSH, inf", "FHO, inf"...
 , 'Location','eastoutside');
xlabel("T^0 behind reflected SW, K");
xticks([2000 3000 4000 5000 6000 7000 8000 9000]);
ylabel("\DeltaT_v/T_v");
title("NO/Ar max \DeltaT_v/T_v");
xlim([2000 9000]);
grid minor
end

%%
if Nerror_plot
%same as Tv_deviation plot
figure
hold on
for var=1:11

j=1;
while (data_experiment(var).n(j,3)~=0)
j=j+1;
end
time_n_exp=data_experiment(var).n(1:j-1,1);
tlim_n=max(time_n_exp);
n_exp=data_experiment(var).n(1:j-1,2);
[time_n_exp, I]=sort(time_n_exp);
n_exp=n_exp(I);

k=1;
for i=1:length(n_exp)
    if time_n_exp(i)<0
        k=i;
    end
end
time_n_exp=time_n_exp(k+1:end);
n_exp=n_exp(k+1:end);

i=2;
while i~=length(time_n_exp)
    if time_n_exp(i)==time_n_exp(i-1)
        time_n_exp=[time_n_exp(1:i-1),
            time_n_exp(i+1:end)];
        n_exp=[n_exp(1:i-1),
            n_exp(i+1:end)];
    end
    i=i+1;
end

i_vibr=1; rel=2;
s1=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH D/6k
s2=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH 3T
s3=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH inf
i_vibr=2;
s4=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO D/6k
s5=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO 3T
s6=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nNO*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO inf

%Uncomment if you want to account cases without relaxation
%dont forgot to fix a legend
% rel=1;
% i_vibr=1;
% s7=spline(data_main(i_vibr,2,var,rel).time, ...
%     data_main(i_vibr,2,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %SSH D/6k
% s8=spline(data_main(i_vibr,3,var,rel).time, ...
%     data_main(i_vibr,3,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %SSH 3T
% s9=spline(data_main(i_vibr,4,var,rel).time, ...
%     data_main(i_vibr,4,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %SSH inf
% i_vibr=2;
% s10=spline(data_main(i_vibr,2,var,rel).time, ...
%     data_main(i_vibr,2,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %FHO D/6k
% s11=spline(data_main(i_vibr,3,var,rel).time, ...
%     data_main(i_vibr,3,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %FHO 3T
% s12=spline(data_main(i_vibr,4,var,rel).time, ...
%     data_main(i_vibr,4,var,rel).nNO*1e3, ...
%     min(time_n_exp):0.1:max(time_n_exp)); %FHO inf


s_exp=makima(time_n_exp, n_exp, min(time_n_exp):0.1:max(time_n_exp)); %EXP

err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
%Uncomment if you want to account cases without relaxation
%dont forgot to fix a legend
% err7=mean(abs(s7-s_exp)./s_exp);
% err8=mean(abs(s8-s_exp)./s_exp);
% err9=mean(abs(s9-s_exp)./s_exp);
% err10=mean(abs(s10-s_exp)./s_exp);
% err11=mean(abs(s11-s_exp)./s_exp);
% err12=mean(abs(s12-s_exp)./s_exp);

p1=plot(TEMP(var), err1, 'ok');
p2=plot(TEMP(var), err4, 'o','color', [0 0.6 0]);
%plot(TEMP(var), err7, 'ob');
%plot(TEMP(var), err10,'o','color', [0.9 0 0]);
p3=plot(TEMP(var), err2, 'sk');
p4=plot(TEMP(var), err5, 's','color',[0 0.6 0]);
%plot(TEMP(var), err8, 'sb');
%plot(TEMP(var), err11, 's','color', [0.9 0 0]);
p5=plot(TEMP(var), err3, 'dk');
p6=plot(TEMP(var), err6,'d','color',[0 0.6 0]);
%plot(TEMP(var), err9, 'db');
%plot(TEMP(var), err12, 'd','color', [0.9 0 0]);
end
[TEMP, I]=sort(TEMP);
exp_err=exp_err(I);
p7=plot(TEMP, exp_err, 'k--');
hold off
set(gca, 'FontName', 'Times New Roman');
legend([p1 p2 p3 p4 p5 p6 p7],"SSH,\it U=D/6k", "FHO,\it U=D/6k"...
 , "SSH,\it U=3T", "FHO,\it U=3T"...
 , "SSH,\it U=\infty", "FHO,\it U=\infty", "mean exp. error"...
 , 'Location','eastoutside');
xlabel("{\it T}^0 behind reflected SW, K");
xticks([2000 3000 4000 5000 6000 7000 8000 9000]);
ylabel("\Delta\it n_{\rm NO}/n_{\rm NO}");
%title("NO/Ar mean \Delta n_{NO}/n_{NO}");
xlim([2000 9000]);
ylim([0 0.8]);
grid minor
end