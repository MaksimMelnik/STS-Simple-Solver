%%Comparison with experimental data and plots
%Streicher's 2022 experiment with NO\O\N\O2\N2\Ar mixtures

clearvars;
%Here you need to load an array with calculated data
load('..\data\Air Streicher experiment\Air_behindRSW_output.mat');
load('..\data\Air Streicher experiment\Air_Streicher_experiment.mat');
info=["21% O2/N2 No.1 P0=64 Torr", "21% O2/N2 No.2 P0=32 Torr"];

Na=6.02214076e23;



%Initialization of calculated data. data_main - the main model, 
% in our case it was a model with included exchange reactions 
% according to the Kunova model. data_minor - is data for comparison 
% with the data_main, in our case these were models without taking into 
% account exchange reactions or with averaging of the Kunova model

data_main=dat1(:,:,:,:,2); %data with i_exch=2 model (full NO spectrum)

%load('..\data\NO_N2 Streicher experiment\NO_N2_behindRSW_output.mat');
data_minor=dat1(:,:,:,:,1); %data with i_exch=1 model(without exch. react.)
data_minor2=dat1(:,:,:,:,3);

%indicators for plots. True if you want to plot a graph,
%False if you dont want to.
%MACRO_plot - plot with macro parameters: T, TvNO, n_NO, n_NO_0
%VDF_NO_plot - plot VDF for NO
%VDF_O2N2_plot - plot VDF for N2 and O2
%TVerror_plot - plot deviation for Tv_NO
%Nerror_plot - plot deviation for n_NO

T_plot=true;
n_plot=true;
VDF_plot=true;
TVerror_plot=false;
Nerror_plot=false;


for var=1:2
%testcases     
%var: %1 - T=4520 K, P=0.393 atm;  2 - 2-N-15 T=6080 K, P=0.097 atm;

i_vibr=2; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on

%Experimental errors and data handling
time_TvO2_err=[];
TvO2_err=[];
err_TvO2=[];
%initialization of experimental data for temperature
j=1;
while (data_experiment(var).TvO2(j,3)~=0)
j=j+1;
end
time_TvO2_exp=data_experiment(var).TvO2(1:j-1,1);
TvO2_exp=data_experiment(var).TvO2(1:j-1,2);
[time_TvO2_exp, I]=sort(time_TvO2_exp);
TvO2_exp=TvO2_exp(I);
tlim_TvO2=max(time_TvO2_exp);
for i=j+1:length(data_experiment(var).TvO2(:,3))-1
    if (data_experiment(var).TvO2(i,3)~=-1)&&...
            (data_experiment(var).TvO2(i-1,3)==data_experiment(var).TvO2(i,3))...
            &&(data_experiment(var).TvO2(i+1,3)==data_experiment(var).TvO2(i,3))

        time_TvO2_err=[time_TvO2_err, (data_experiment(var).TvO2(i,1)+...
           data_experiment(var).TvO2(i-1,1)+data_experiment(var).TvO2(i+1,1))/3];
        TvO2_err=[TvO2_err, data_experiment(var).TvO2(i-1,2)];
        err_TvO2=[err_TvO2, (data_experiment(var).TvO2(i,2) -...
            data_experiment(var).TvO2(i+1,2))/2];
    end
end

if var==2
time_TvN2_exp=data_experiment(var).TvN2(:,1);
TvN2_exp=data_experiment(var).TvN2(:,2);
[time_TvN2_exp, I]=sort(time_TvN2_exp);
TvN2_exp=TvN2_exp(I);

time_TvN2_dsmc=data_experiment(var).TvN2_dsmc(:,1);
TvN2_dsmc=data_experiment(var).TvN2_dsmc(:,2);
[time_TvN2_dsmc, I]=sort(time_TvN2_dsmc);
TvN2_dsmc=TvN2_dsmc(I);
end

time_TvO2_dsmc=data_experiment(var).TvO2_dsmc(:,1);
TvO2_dsmc=data_experiment(var).TvO2_dsmc(:,2);
[time_TvO2_dsmc, I]=sort(time_TvO2_dsmc);
TvO2_dsmc=TvO2_dsmc(I);


time_T_exp=data_experiment(var).T(:,1);
T_exp=data_experiment(var).T(:,2);
[time_T_exp, I]=sort(time_T_exp);
T_exp=T_exp(I);

time_T_dsmc=data_experiment(var).T_dsmc(:,1);
T_dsmc=data_experiment(var).T_dsmc(:,2);
[time_T_dsmc, I]=sort(time_T_dsmc);
T_dsmc=T_dsmc(I);
%End of experimental data handling

%Macro parameters plot

if T_plot
figure("Position", [0, 0, 900, 800])
t=tiledlayout(2, 2, "TileSpacing", "compact");
title(t, "Case " + info(var), 'FontName',...
 'Times New Roman');

nexttile %Oxygen vibr temperature TvO2
set(gca, 'FontName', 'Times New Roman');
hold on
errorbar(time_TvO2_err, TvO2_err, err_TvO2, 'sqk', 'MarkerFaceColor',...
    'k','MarkerSize',1, 'LineWidth', 1)
time_TvO2_err=[0 time_TvO2_err tlim_TvO2];
ppp=spline(time_TvO2_exp, TvO2_exp, time_TvO2_err);
p4=plot(time_TvO2_err, ppp, 'k-', 'LineWidth', 2, 'DisplayName', ...
    "\it T -\rm experiment ");
p5=plot(time_TvO2_dsmc, TvO2_dsmc, 'sk', 'LineWidth', 1, 'DisplayName', ...
    "\it T -\rm DSMC");
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).TvO2,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k with exch." );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).TvO2,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T with exch." );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).TvO2,...
    'm-','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty with exch." );

legend([p1 p2 p3 p4 p5],'Location','best');
xlim([-5 60]);
ylim([0 max(TvO2_exp)+1000]);
xlabel("\it t,\rm \mus");
ylabel("\it {\it T}_v^{O_2},\rm K");
hold off
grid minor

nexttile %Vibrational temperature of N2
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).TvN2,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it T}_v^{N_2} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).TvN2,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it T}_v^{N_2} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).TvN2,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it T}_v^{N_2} -\it U=\infty " );

if var==2
p4=plot(time_TvN2_exp, TvN2_exp, 'k-', 'LineWidth', 2, 'DisplayName',...
    "{\it T}_v^{N_2} - experiment");
p5=plot(time_TvN2_dsmc, TvN2_dsmc, 'sk', 'LineWidth', 1, 'DisplayName',...
    "{\it T}_v^{N_2} - DSMC");
end
xlim([-5 60]);
%ylim([0 6000]);
if var==2
legend([p1 p2 p3 p4 p5],'Location','best');
else
legend([p1 p2 p3],'Location','best');
end
xlabel("{\it t}, \mus");
ylabel("{\it T}_v^{N_2}, K");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).TvNO,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it T}_v^{NO} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).TvNO,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it T}_v^{NO} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).TvNO,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it T}_v^{NO} -\it U=\infty " );


xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("{\it T}_v^{NO}, K");
hold off
grid minor


nexttile %Translational temp
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).T,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it T} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).T,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it T} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).T,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it T} -\it U=\infty " );
p4=plot(time_T_exp, T_exp, 'k-', 'LineWidth', 2, 'DisplayName',...
    "{\it T} - experiment");
p5=plot(time_T_dsmc, T_dsmc, 'sk', 'LineWidth', 1, 'DisplayName',...
    "{\it T} - DSMC");
xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3 p4 p5],'Location','best');
xlabel("{\it t}, \mus");
ylabel("{\it T}, K");
hold off
grid minor


%exportgraphics(t, "NON2Ar_Scanlon_case"+num2str(var)+".jpg");
end

if n_plot
figure("Position", [0, 0, 1500, 800])
t=tiledlayout(2, 3, "TileSpacing", "compact");
title(t, "Case " + info(var), 'FontName',...
 'Times New Roman');

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nNO*1e3,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{NO} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nNO*1e3,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{NO} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nNO*1e3,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{NO} -\it U=\infty " );

p11=plot(data_minor(i_vibr,2,var,rel).time, data_minor(i_vibr,2,var,rel).nNO*1e3,...
    'r--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k without exch." );
p22=plot(data_minor(i_vibr,3,var,rel).time, data_minor(i_vibr,3,var,rel).nNO*1e3,...
    'b--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T without exch." );
p33=plot(data_minor(i_vibr,4,var,rel).time, data_minor(i_vibr,4,var,rel).nNO*1e3,...
    'm--','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty without exch." );

p111=plot(data_minor2(i_vibr,2,var,rel).time, data_minor2(i_vibr,2,var,rel).nNO*1e3,...
    'r-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k with exch. NO(0)" );
p222=plot(data_minor2(i_vibr,3,var,rel).time, data_minor2(i_vibr,3,var,rel).nNO*1e3,...
    'b-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T with exch. NO(0)" );
p333=plot(data_minor2(i_vibr,4,var,rel).time, data_minor2(i_vibr,4,var,rel).nNO*1e3,...
    'm-.','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty with exch. NO(0)" );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("n_{NO}, mmol/m^3");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
% p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nN2./(data_main(i_vibr,2,var,rel).nNO+data_main(i_vibr,2,var,rel).nN2 +data_main(i_vibr,2,var,rel).nO2 + data_main(i_vibr,2,var,rel).nN + data_main(i_vibr,2,var,rel).nO),...
%     'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=D/6k " );
% p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nN2./(data_main(i_vibr,3,var,rel).nNO + data_main(i_vibr,3,var,rel).nN2 + data_main(i_vibr,3,var,rel).nO2 + data_main(i_vibr,3,var,rel).nN + data_main(i_vibr,3,var,rel).nO),...
%     'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=3T " );
% p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nN2./(data_main(i_vibr,4,var,rel).nNO + data_main(i_vibr,4,var,rel).nN2 + data_main(i_vibr,4,var,rel).nO2 + data_main(i_vibr,4,var,rel).nN + data_main(i_vibr,4,var,rel).nO),...
%     'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=\infty " );
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nN2*1e3,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nN2*1e3,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nN2*1e3,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{N_2} -\it U=\infty " );

p11=plot(data_minor(i_vibr,2,var,rel).time, data_minor(i_vibr,2,var,rel).nN2*1e3,...
    'r--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k without exch." );
p22=plot(data_minor(i_vibr,3,var,rel).time, data_minor(i_vibr,3,var,rel).nN2*1e3,...
    'b--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T without exch." );
p33=plot(data_minor(i_vibr,4,var,rel).time, data_minor(i_vibr,4,var,rel).nN2*1e3,...
    'm--','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty without exch." );

p111=plot(data_minor2(i_vibr,2,var,rel).time, data_minor2(i_vibr,2,var,rel).nN2*1e3,...
    'r-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k with exch. NO(0)" );
p222=plot(data_minor2(i_vibr,3,var,rel).time, data_minor2(i_vibr,3,var,rel).nN2*1e3,...
    'b-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T with exch. NO(0)" );
p333=plot(data_minor2(i_vibr,4,var,rel).time, data_minor2(i_vibr,4,var,rel).nN2*1e3,...
    'm-.','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty with exch. NO(0)" );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("n_{NO}, mmol/m^3");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).p,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it p} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).p,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it p} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).p,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it p} -\it U=\infty " );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("\it p, \rm Torr");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
% p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nO2./(data_main(i_vibr,2,var,rel).nNO+data_main(i_vibr,2,var,rel).nN2 +data_main(i_vibr,2,var,rel).nO2 + data_main(i_vibr,2,var,rel).nN + data_main(i_vibr,2,var,rel).nO),...
%     'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=D/6k " );
% p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nO2./(data_main(i_vibr,3,var,rel).nNO + data_main(i_vibr,3,var,rel).nN2 + data_main(i_vibr,3,var,rel).nO2 + data_main(i_vibr,3,var,rel).nN + data_main(i_vibr,3,var,rel).nO),...
%     'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=3T " );
% p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nO2./(data_main(i_vibr,4,var,rel).nNO + data_main(i_vibr,4,var,rel).nN2 + data_main(i_vibr,4,var,rel).nO2 + data_main(i_vibr,4,var,rel).nN + data_main(i_vibr,4,var,rel).nO),...
%     'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=\infty " );

p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nO2*1e3,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nO2*1e3,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nO2*1e3,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{O_2} -\it U=\infty " );

p11=plot(data_minor(i_vibr,2,var,rel).time, data_minor(i_vibr,2,var,rel).nO2*1e3,...
    'r--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k without exch." );
p22=plot(data_minor(i_vibr,3,var,rel).time, data_minor(i_vibr,3,var,rel).nO2*1e3,...
    'b--', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T without exch." );
p33=plot(data_minor(i_vibr,4,var,rel).time, data_minor(i_vibr,4,var,rel).nO2*1e3,...
    'm--','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty without exch." );

p111=plot(data_minor2(i_vibr,2,var,rel).time, data_minor2(i_vibr,2,var,rel).nO2*1e3,...
    'r-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=D/6k with exch. NO(0)" );
p222=plot(data_minor2(i_vibr,3,var,rel).time, data_minor2(i_vibr,3,var,rel).nO2*1e3,...
    'b-.', 'LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=3T with exch. NO(0)" );
p333=plot(data_minor2(i_vibr,4,var,rel).time, data_minor2(i_vibr,4,var,rel).nO2*1e3,...
    'm-.','LineWidth', 1.5, 'DisplayName', "\it {\it T}_v^{O_2} - U=\infty with exch. NO(0)" );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("n_{NO}, mmol/m^3");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
% p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nN./(data_main(i_vibr,2,var,rel).nNO+data_main(i_vibr,2,var,rel).nN2 +data_main(i_vibr,2,var,rel).nO2 + data_main(i_vibr,2,var,rel).nN + data_main(i_vibr,2,var,rel).nO),...
%     'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=D/6k " );
% p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nN./(data_main(i_vibr,3,var,rel).nNO + data_main(i_vibr,3,var,rel).nN2 + data_main(i_vibr,3,var,rel).nO2 + data_main(i_vibr,3,var,rel).nN + data_main(i_vibr,3,var,rel).nO),...
%     'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=3T " );
% p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nN./(data_main(i_vibr,4,var,rel).nNO + data_main(i_vibr,4,var,rel).nN2 + data_main(i_vibr,4,var,rel).nO2 + data_main(i_vibr,4,var,rel).nN + data_main(i_vibr,4,var,rel).nO),...
%     'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=\infty " );

p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nN,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nN,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nN,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{N} -\it U=\infty " );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("n_{N}, mmol/m^3");
hold off
grid minor

nexttile %Vibrational temperature of NO
set(gca, 'FontName', 'Times New Roman');
hold on
% p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nO./(data_main(i_vibr,2,var,rel).nNO+data_main(i_vibr,2,var,rel).nN2 +data_main(i_vibr,2,var,rel).nO2 + data_main(i_vibr,2,var,rel).nN + data_main(i_vibr,2,var,rel).nO),...
%     'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=D/6k " );
% p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nO./(data_main(i_vibr,3,var,rel).nNO + data_main(i_vibr,3,var,rel).nN2 + data_main(i_vibr,3,var,rel).nO2 + data_main(i_vibr,3,var,rel).nN + data_main(i_vibr,3,var,rel).nO),...
%     'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=3T " );
% p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nO./(data_main(i_vibr,4,var,rel).nNO + data_main(i_vibr,4,var,rel).nN2 + data_main(i_vibr,4,var,rel).nO2 + data_main(i_vibr,4,var,rel).nN + data_main(i_vibr,4,var,rel).nO),...
%     'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=\infty " );

p1=plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).nO,...
    'r-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=D/6k " );
p2=plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).nO,...
    'b-', 'LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=3T " );
p3=plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).nO,...
    'm-','LineWidth', 1.5, 'DisplayName', "{\it n}_{O} -\it U=\infty " );

xlim([-5 60]);
%ylim([0 6000]);
legend([p1 p2 p3],'Location','best');
xlabel("{\it t}, \mus");
ylabel("n_{O}, mmol/m^3");
hold off
grid minor
end

%% VDF plots
if VDF_plot
addpath('../src/');
addpath('../data/')
load('particles.mat', "NO", "N", "O", "Ar", "O2", "N2");
k=1.380649e-23; Na=6.02214076e23;
NO.num_elex_levels=1;
N2.num_elex_levels=1;
O2.num_elex_levels=1;
rel=2;
i_vibr=2;

i_U=2; %Choose dissociation model: 2-D/6k, 3-3T, 4-infty

[~,j1]=min(abs(data_main(i_vibr, i_U, var, rel).time - tlim_TvO2/3));
%determination of the moment in time tlim/3
[~,j2]=min(abs(data_main(i_vibr, i_U, var, rel).time - 2*tlim_TvO2/3));
%determination of the moment in time 2*tlim/3
[~,j3]=min(abs(data_main(i_vibr, i_U, var, rel).time - tlim_TvO2));
%determination of the moment in time tlim
Na=6.02214076e23; %Avogadro constant

J=[5 j1 j2 j3];
%total number density of the mixture in case without exch. reactions
nSUMminor=(data_minor(i_vibr, i_U, var, rel).nNO + ...
    data_minor(i_vibr, i_U, var, rel).nN2 + ...
    data_minor(i_vibr, i_U, var, rel).nO2 + ...
    data_minor(i_vibr, i_U, var, rel).nN + ...
    data_minor(i_vibr, i_U, var, rel).nO)*Na;

%total number density of the mixture in case with exch. reactions
nSUMmain=(data_main(i_vibr, i_U, var, rel).nNO + ...
    data_main(i_vibr, i_U, var, rel).nN2 + ...
    data_main(i_vibr, i_U, var, rel).nO2 + ...
    data_main(i_vibr, i_U, var, rel).nN + ...
    data_main(i_vibr, i_U, var, rel).nO)*Na;

[~,j11]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - tlim_TvO2/3));
[~,j22]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - 2*tlim_TvO2/3));
[~,j33]=min(abs(data_minor(i_vibr, i_U, var, rel).time...
    - tlim_TvO2));

[~,j111]=min(abs(data_minor2(i_vibr, i_U, var, rel).time...
    - tlim_TvO2/3));
[~,j222]=min(abs(data_minor2(i_vibr, i_U, var, rel).time...
    - 2*tlim_TvO2/3));
[~,j333]=min(abs(data_minor2(i_vibr, i_U, var, rel).time...
    - tlim_TvO2));
JJ=[5 j11 j22 j33];
JJJ=[5 j111 j222 j333];
tlim = [0 tlim_TvO2/3 2*tlim_TvO2/3 tlim_TvO2];

for i=1:4
n_bolzNO=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(J(i)),...
    data_main(i_vibr, i_U, var, rel).nNO(J(i))*Na, NO);
% n_bolzNO_2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j1),...
%     data_main(i_vibr, i_U, var, rel).nNO(j1)*Na, NO);
% n_bolzNO_3=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j2),...
%     data_main(i_vibr, i_U, var, rel).nNO(j2)*Na, NO);
% n_bolzNO_4=density_f_exc(data_main(i_vibr, i_U, var, rel).TvNO(j3),...
%     data_main(i_vibr, i_U, var, rel).nNO(j3)*Na, NO);

n_bolzO2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvO2(J(i)),...
    data_main(i_vibr, i_U, var, rel).nO2(J(i))*Na, O2);
% n_bolzO2_2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvO2(j1),...
%     data_main(i_vibr, i_U, var, rel).nO2(j1)*Na, O2);
% n_bolzO2_3=density_f_exc(data_main(i_vibr, i_U, var, rel).TvO2(j2),...
%     data_main(i_vibr, i_U, var, rel).nO2(j2)*Na, O2);
% n_bolzO2_4=density_f_exc(data_main(i_vibr, i_U, var, rel).TvO2(j3),...
%     data_main(i_vibr, i_U, var, rel).nO2(j3)*Na, O2);

n_bolzN2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvN2(J(i)),...
    data_main(i_vibr, i_U, var, rel).nN2(J(i))*Na, N2);
% n_bolzN2_2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvN2(j1),...
%     data_main(i_vibr, i_U, var, rel).nN2(j1)*Na, N2);
% n_bolzN2_3=density_f_exc(data_main(i_vibr, i_U, var, rel).TvN2(j2),...
%     data_main(i_vibr, i_U, var, rel).nN2(j2)*Na, N2);
% n_bolzN2_4=density_f_exc(data_main(i_vibr, i_U, var, rel).TvN2(j3),...
%     data_main(i_vibr, i_U, var, rel).nN2(j3)*Na, N2);

figure("Position", [200, 0, 1300, 500])

% ylim_low=min([min(n_bolz1./nSUMmain(1)),...
%     min(n_bolz2./nSUMmain(j1)),...
%     min(n_bolz3./nSUMmain(j2)),...
%     min(n_bolz4./nSUMmain(j3)),...
%     min(data_main(i_vibr, i_U, var, rel).ni_NO(1, :)./...
%     nSUMmain(1)), min(data_main(i_vibr, i_U, var, rel).ni_NO(j1, :)./...
%     nSUMmain(j1)), min(data_main(i_vibr, i_U, var, rel).ni_NO(j2, :)./...
%     nSUMmain(j2)), min(data_main(i_vibr, i_U, var, rel).ni_NO(j3, :)./...
%     nSUMmain(j3)), min(data_minor(i_vibr, i_U, var, rel).ni_NO(1, :)./...
%     nSUMminor(1)), min(data_minor(i_vibr, i_U, var, rel).ni_NO(j11, :)./...
%     nSUMminor(j11)), min(data_minor(i_vibr, i_U, var, rel).ni_NO(j22, :)./...
%     nSUMminor(j22)), min(data_minor(i_vibr, i_U, var, rel).ni_NO(j33, :)./...
%     nSUMminor(j33))])*1e-1;
ylim_low=1e-20;
set(gca, 'FontName', 'Times New Roman');
t=tiledlayout(1, 3, "TileSpacing", "compact");
%title(t, "Case " + info(var) + ". T^0=" + num2str(round(max((T_exp))))...
% + " K", 'FontName', 'Times New Roman');

nexttile 
semilogy(0:NO.num_vibr_levels(1)-1, n_bolzNO./nSUMmain(J(i)), 'k-',...
    'LineWidth',2.0, "DisplayName","Boltzmann");
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:NO.num_vibr_levels(1)-1,...
    data_minor(i_vibr, i_U, var, rel).ni_NO(JJ(i), :)./...
    nSUMminor(J(i)) ,'-','color',[0 0.6 0],'LineWidth',2.0,...
    "DisplayName","Without exch. react.");
semilogy(0:NO.num_vibr_levels(1)-1, ...
    data_main(i_vibr, i_U, var, rel).ni_NO(J(i), :)./...
    nSUMmain(J(i)) , 'm-','LineWidth',2.0, "DisplayName",...
    "With exch. react.");
% semilogy(0:NO.num_vibr_levels(1)-1, ...
%     data_minor2(i_vibr, i_U, var, rel).ni_NO(JJJ(i), :)./...
%     nSUMmain(J(i)) , 'g-','LineWidth',2.0, "DisplayName",...
%     "With exch. react. NO(0)");

xlabel('NO vibr. level,\it i');
ylabel('\it n_{\rm NO\iti}/\itn');
%lgd=legend('Location','best');
xlim([0 NO.num_vibr_levels(1)-1]);
ylim([-inf 1])
text(10,1e-1,"{\it t} = " + num2str(round(tlim(i))) + "\mus," + " {\it T}_v^{NO}({\itt}) = "...
    +num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(J(i))))+" K",...
    'FontSize',12, 'FontName', 'Times New Roman');
hold off
box off
grid minor

nexttile %t=tlim/3
semilogy(0:O2.num_vibr_levels(1)-1, n_bolzO2./nSUMmain(J(i)), 'k-', ...
    'LineWidth',2.0, "DisplayName","Boltzmann");
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:O2.num_vibr_levels(1)-1,...
    data_minor(i_vibr, 2, var, rel).ni_O2(JJ(i), :)./...
    nSUMminor(J(i)) ,'-','color',[0 0.6 0],'LineWidth',2.0, "DisplayName","Without exch. react.");
semilogy(0:O2.num_vibr_levels(1)-1,...
    data_main(i_vibr, 2, var, rel).ni_O2(J(i), :)./...
    nSUMmain(J(i)) , 'm-','LineWidth',2.0, "DisplayName","With exch. react.");
semilogy(0:O2.num_vibr_levels(1)-1, ...
    data_minor2(i_vibr, i_U, var, rel).ni_O2(JJJ(i), :)./...
    nSUMmain(J(i)) , '-','color', [0.9 0 0], 'LineWidth',2.0, "DisplayName","With exch. react. NO(0)");

xlabel('O_2 vibr. level,\it i');
ylabel('\it n_{\rm O_2\iti}/\itn');
xlim([0 O2.num_vibr_levels(1)-1]);
%ylim([ylim_low 1])
text(10,1e-1,"{\it t} = "+ num2str(round(tlim(i))) + " \mus,"  + ...
    " {\it T}_v^{O_2}({\itt}) = "+ ...
    num2str(round(data_main(i_vibr, i_U, var, rel).TvO2(J(i))))+" K", ...
    'FontSize',12, 'FontName', 'Times New Roman')
lgd=legend('Location','south');
hold off
box off
grid minor

nexttile %t=2*tlim/3
semilogy(0:N2.num_vibr_levels(1)-1, n_bolzN2./nSUMmain(J(i)), ...
    'k-', 'LineWidth',2.0);
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times New Roman');
hold on;
semilogy(0:N2.num_vibr_levels(1)-1, ...
    data_minor(i_vibr, i_U, var, rel).ni_N2(JJ(i), :)./ ...
    nSUMminor(J(i)) ,'-','color',[0 0.6 0],'LineWidth',2.0);
semilogy(0:N2.num_vibr_levels(1)-1, ...
    data_main(i_vibr, i_U, var, rel).ni_N2(J(i), :)./ ...
    nSUMmain(J(i)) , 'm-','LineWidth',2.0);
semilogy(0:N2.num_vibr_levels(1)-1, ...
    data_minor2(i_vibr, i_U, var, rel).ni_N2(JJJ(i), :)./...
    nSUMmain(J(i)) , '-','color', [0.9 0 0],'LineWidth',2.0);
xlabel('NO vibr. level,\it i');
ylabel('\it n_{\rm NO\iti}/\itn');
xlim([0 N2.num_vibr_levels(1)-1]);
%ylim([ylim_low 1])
text(10,1e-1,"{\it t} = "+ num2str(round(tlim(i))) + " \mus," ...
    + " {\it T}_v^{N_2}({\itt}) = "+ ...
    num2str(round(data_main(i_vibr, i_U, var, rel).TvN2(J(i))))+" K", ...
    'FontSize',12, 'FontName', 'Times New Roman')
hold off
box off
grid minor

lgd.Layout.Tile = 'south';
end
% nexttile %t=tlim
% semilogy(0:NO.num_vibr_levels(1)-1, n_bolz4./nSUMmain(j3), 'k-', ...
%     'LineWidth',2.0);
% set(gca,'FontSize',12);
% set(gca, 'FontName', 'Times New Roman');
% hold on;
% semilogy(0:NO.num_vibr_levels(1)-1, ...
%     data_minor(i_vibr, i_U, var, rel).ni_NO(j33, :)./ ...
%     nSUMminor(j33) ,'-','color',[0 0.6 0],'LineWidth',2.0);
% semilogy(0:NO.num_vibr_levels(1)-1, ...
%     data_main(i_vibr, i_U, var, rel).ni_NO(j3, :)./ ...
%     nSUMmain(j3) , 'm-','LineWidth',2.0);
% xlabel('NO vibr. level,\it i', 'FontName','Times New Roman');
% ylabel('\it n_{\rm NO\iti}/\itn');
% xlim([0 NO.num_vibr_levels(1)-1]);
% ylim([ylim_low 1])
% text(3, ylim_low*1e5,"{\it t} = "+num2str(round(tlim_TvO2))+" \mus," ...
%     + " {\it T}_v^{NO}({\itt}) = "+ ...
%     num2str(round(data_main(i_vibr, i_U, var, rel).TvNO(j3)))+" K", ...
%     'FontSize',12, 'FontName', 'Times New Roman')
% hold off
% box off
% grid minor
% lgd.Layout.Tile = 'south';


%exportgraphics(t, "NON2Ar_VDF_case"+num2str(var)+".jpg");
end



% if VDF_plot
% figure("Position", [200, 0, 900, 500])
% n_bolz_N2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvN2(j33), ...
%     data_main(i_vibr, i_U, var, rel).nN2(j33)*Na, N2);
% n_bolz_O2=density_f_exc(data_main(i_vibr, i_U, var, rel).TvO2(j33), ...
%     data_main(i_vibr, i_U, var, rel).nO2(j33)*Na, O2);
% 
% ylim_low=min([min(n_bolz_N2./nSUMmain(j33)),...
% min(n_bolz_O2./nSUMmain(j33)), ...
% min(data_minor(i_vibr, i_U, var, rel).ni_N2(j33, :)./nSUMminor(j33)), ...
% min(data_main(i_vibr, i_U, var, rel).ni_N2(j3, :)./nSUMmain(j3)),...
% min(data_minor(i_vibr, i_U, var, rel).ni_O2(j33, :)./nSUMminor(j33)), ...
% min(data_main(i_vibr, i_U, var, rel).ni_O2(j3, :)./nSUMmain(j3))])*1e-1;
% 
% set(gca, 'FontName', 'Times New Roman');
% t=tiledlayout(1, 2, "TileSpacing", "compact");
% %title(t, "Case " + info(var) + ". T^0=" + num2str(round(max((T_exp))))+...
% % " K", 'FontName', 'Times New Roman');
% 
% nexttile
% semilogy(0:N2.num_vibr_levels(1)-1, n_bolz_N2./nSUMmain(j33), ...
%     'k-', 'LineWidth',2.0, "DisplayName","Boltzmann");
% set(gca,'FontSize',12);
% set(gca, 'FontName', 'Times New Roman');
% hold on;
% semilogy(0:N2.num_vibr_levels(1)-1, ...
%     data_minor(i_vibr, i_U, var, rel).ni_N2(j33, :)./ ...
%     nSUMminor(j33),'-','color',[0 0.6 0],'LineWidth',2.0, ...
%     "DisplayName","Without exch. react.");
% semilogy(0:N2.num_vibr_levels(1)-1, ...
%     data_main(i_vibr, i_U, var, rel).ni_N2(j3, :)./ ...
%     nSUMmain(j3) , 'm-','LineWidth',2.0, ...
%     "DisplayName","With exch. react. and vibr. activation of react. prod. (NO)");
% xlabel('N_2 vibr. level,\it i', 'FontName','Times New Roman');
% ylabel('\it n_{\rm N_2\iti}/\itn');
% xlim([0 N2.num_vibr_levels(1)-1]);
% ylim([ylim_low 1]);
% text(3, ylim_low*1e1,"{\it t} = "+num2str(round(tlim_TvO2))+" \mus," + ...
%     " {\it T}_v^{N_2}({\itt}) = "+ ...
%     num2str(round(data_main(i_vibr, i_U, var, rel).TvN2(j3)))+" K", ...
%     'FontSize',12, 'FontName', 'Times New Roman')
% lgd=legend('Location','south');
% hold off
% box off
% grid minor
% 
% nexttile
% semilogy(0:O2.num_vibr_levels(1)-1, n_bolz_O2./nSUMmain(j33), ...
%     'k-', 'LineWidth',2.0);
% set(gca,'FontSize',12);
% set(gca, 'FontName', 'Times New Roman');
% hold on;
% semilogy(0:O2.num_vibr_levels(1)-1, ...
%     data_minor(i_vibr, i_U, var, rel).ni_O2(j33, :)./ ...
%     nSUMminor(j33) ,'-','color',[0 0.6 0],'LineWidth',2.0);
% semilogy(0:O2.num_vibr_levels(1)-1, ...
%     data_main(i_vibr, i_U, var, rel).ni_O2(j3, :)./ ...
%     nSUMmain(j33) , 'm-','LineWidth',2.0);
% xlabel('O_2 vibr. level,\it i', 'FontName','Times New Roman');
% ylabel('\it n_{\rm O_2\iti}/\itn');
% xlim([0 O2.num_vibr_levels(1)-1]);
% ylim([ylim_low 1]);
% text(3, ylim_low*1e1,"{\it t} = "+num2str(round(tlim_TvO2))+" \mus," + ...
%     " {\it T}_v^{O_2}({\itt}) = "+ ...
%     num2str(round(data_main(i_vibr, i_U, var, rel).TvO2(j3)))+" K", ...
%     'FontSize',12, 'FontName', 'Times New Roman')
% hold off
% box off
% grid minor
% lgd.Layout.Tile = 'south';
% end
% 
% rmpath('../src/')
% rmpath('../data/')
 end


%%
TEMP=[3580 6080 6970 1880 3570 5190 3730 6240 8210 1920 3490 6050];
exp_err=[0.022 0.30 0.25 0.01 0.025 0.28 0.014 0.20 0.48 0.015 0.03 0.3];
if Nerror_plot
figure
hold on
for var=1:12
  
%separation of data and errors
j=1;
while (data_experiment(var).n(j,3)~=0)
j=j+1;
end
time_n_exp=data_experiment(var).n(1:j-1,1);
tlim_TvO2=max(time_n_exp);
n_exp=data_experiment(var).n(1:j-1,2);
[time_n_exp, I]=sort(time_n_exp);
n_exp=n_exp(I);

%throwing out points with negative time
k=1;
for i=1:length(n_exp)
    if time_n_exp(i)<0
        k=i;
    end
end
time_n_exp=time_n_exp(k+1:end);
n_exp=n_exp(k+1:end);

%throwing out points with the same abscissa
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

%interpolation of calculated data for experimental time
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

%calculating of mean deviation
err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
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
xticks([1000 2000 3000 4000 5000 6000 7000 8000 9000]);
ylabel("\Delta \it n_{\rm NO}/n_{\rm NO}");
%title("NO/N2 mean \Delta n_{NO}/n_{NO}");
xlim([1000 9000]);
grid minor
end


%%
if TVerror_plot

figure
hold on
for var=1:12

time_TvNO_exp=data_experiment(var).TvNO(:,1);
TvNO_exp=data_experiment(var).TvNO(:,2);
[time_TvNO_exp, I]=sort(time_TvNO_exp);
TvNO_exp=TvNO_exp(I);
k=1;
for i=1:length(TvNO_exp)
    if time_TvNO_exp(i)<0
        k=i;
    end
end
time_TvNO_exp=time_TvNO_exp(k+1:end);
TvNO_exp=TvNO_exp(k+1:end);

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

err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
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
xticks([1000 2000 3000 4000 5000 6000 7000 8000 9000]);
ylabel("\DeltaT_v/T_v");
title("NO/N2 max \DeltaT_v/T_v");
xlim([1000 9000]);
grid minor
end