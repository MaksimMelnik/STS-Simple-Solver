%%Comparison with experimental data and plots
%Streicher's 2021 experiment with O2\O\Ar mixtures
clearvars;
%Here you need to load an array with calculated data
load('..\data\O2_Ar Streicher experiment\O2Ar_behindRSW_output.mat');
load('..\data\O2_Ar Streicher experiment\O2_Ar_Streicher21_experiment.mat');
info=["50% No.1 (03)", "No.2: O_2 - 50%; Ar - 50%", "50% No.3 (14)" ,...
    "20% No.1 (02)" ,"20% No.2 (08)", "20% No.3 (14)",...
    "100% No.1 (01)","100% No.2 (06)","100% No.3 (08)"];

%Initialization of calculated data. data_main - the main model, 
% in our case it was a model with included VV,VT, diss-rec processes 
data_main=dat1;
% here you need to load a minor data, if you want to compare new model with 
% main data
data_minor=dat1;

%indicators for plots. True if you want to plot a graph,
%False if you dont want to.
%MACRO_plot - plot with macro parameters: TvO2, n_O2
%P_plot - plot with pressure (only for testcase#1)
%TVerror_plot - plot deviation for Tv_NO
%Nerror_plot - plot deviation for n_NO

MACRO_plot=true;
P_plot=true;
TVerror_plot=false;
Nerror_plot=false;

for var=1%:9
%testcases
%var: %1 - 50-03 T=8110 K, P=75 Torr;  2 - 50-11 T=10470 K, P=53 Torr; 
%3 - 50-13 T=11410 K, P=30 Torr; 4 - 20-02 T=7840 K, P=130 Torr;
%5 - 20-08 T=10310 K, P=97 Torr; 6 - 20-14 T=13830 K, P=33 Torr;
%7 - 100-01 T=6230 K, P=57 Torr; 8 - 100-06 T=7940 K, P=41 Torr;
%9 - 100-08 T=9560 K, P=34 Torr;

i_vibr=1; %model of vibrational enegry exchange 1 - SSH, 2 - FHO
rel=2; %switcher of relaxation between SWs: 1 - off, 2 - on


%Experimental errors and data handling
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
    if (data_experiment(var).Tv(i,3)~=-1)&&...
        (data_experiment(var).Tv(i-1,3)==data_experiment(var).Tv(i,3))&&...
        (data_experiment(var).Tv(i+1,3)==data_experiment(var).Tv(i,3))
        time_Tv_err=[time_Tv_err, (data_experiment(var).Tv(i,1)+ ...
            data_experiment(var).Tv(i-1,1)+ ...
            data_experiment(var).Tv(i+1,1))/3];
        Tv_err=[Tv_err, data_experiment(var).Tv(i-1,2)];
        err_Tv=[err_Tv, (data_experiment(var).Tv(i,2) - ...
            data_experiment(var).Tv(i+1,2))/2];
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
%End of experimental data handling



if MACRO_plot
%%Temperature & Number density
figure("Position", [0, 0, 900, 400])
t=tiledlayout(1, 2, "TileSpacing", "compact");
%title(t, "Case " + info(var)+"; vibr. model: SSH", 'FontName',...
% 'Times New Roman');

%Uncomment if you want to add translational temperature
%dont forgot to fix tiledlayout!!

% nexttile
% hold on
% plot(time_T_exp, T_exp, 'k--', 'LineWidth', 1.5, 'DisplayName', ...
%     "T - exp case " + num2str(var));
% plot(data_main(i_vibr,2,var,rel).time, data_main(i_vibr,2,var,rel).T, ...
%     'r--', 'LineWidth', 1.5, 'DisplayName', "T - U=D/6k " );
% plot(data_main(i_vibr,3,var,rel).time, data_main(i_vibr,3,var,rel).T, ...
%     'b--', 'LineWidth', 1.5, 'DisplayName', "T - U=3T " );
% plot(data_main(i_vibr,4,var,rel).time, data_main(i_vibr,4,var,rel).T, ...
%     'm--','LineWidth', 1.5, 'DisplayName', "T - U=inf " );
% legend('Location','best');
% xlim([0 50]);
% ylim([0 14000]);
% xlabel("t, \mus");
% ylabel("T, K");
% hold off
% grid minor

nexttile
hold on
set(gca, 'FontName', 'Times New Roman');
p1=plot(time_Tv_exp, Tv_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', ...
    "\it T_{\rm v}^{\rm O_2}\rm - experiment ");
p2=plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).Tv,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=D/6k " );
p3=plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).Tv,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=3T " );
p4=plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).Tv, 'm-','LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=\infty " );
p5=plot(data_minor(i_vibr,2,var,rel).time, ...
    data_minor(i_vibr,2,var,rel).Tv,'r--', 'LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=D/6k minor" );
p6=plot(data_minor(i_vibr,3,var,rel).time, ...
    data_minor(i_vibr,3,var,rel).Tv,'b--', 'LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=3T minor" );
p7=plot(data_minor(i_vibr,4,var,rel).time, ...
    data_minor(i_vibr,4,var,rel).Tv, 'm--','LineWidth', ...
    1.5, 'DisplayName', "\it T_{\rm v}^{\rm O_2} - U=\infty minor" );
errorbar(time_Tv_err, Tv_err, err_Tv, 'sqk', ...
    'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
xlim([0 max(time_Tv_exp)]);
legend([p1 p2 p3 p4 p5 p6 p7],'Location','best');
xlabel("\it t\rm, \mus");
ylabel("\it T_{\rm v}^{\rm O_2}\rm, K");
hold off
grid minor


nexttile 
hold on
set(gca, 'FontName', 'Times New Roman');
p1=plot(time_n_exp, n_exp, 'k-', 'LineWidth', 1.5, 'DisplayName', ...
    "\it n_{\rm O_2}\rm - experiment ");
p2=plot(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nO2*1e3,'r-', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=D/6k " );
p3=plot(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nO2*1e3,'b-', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=3T " );
p4=plot(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nO2*1e3, 'm-','LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=\infty " );
p5=plot(data_minor(i_vibr,2,var,rel).time, ...
    data_minor(i_vibr,2,var,rel).nO2*1e3,'r--', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=D/6k minor" );
p6=plot(data_minor(i_vibr,3,var,rel).time, ...
    data_minor(i_vibr,3,var,rel).nO2*1e3,'b--', 'LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=3T minor" );
p7=plot(data_minor(i_vibr,4,var,rel).time, ...
    data_minor(i_vibr,4,var,rel).nO2*1e3, 'm--','LineWidth', ...
    1.5, 'DisplayName', "\it n_{\rm O_2} - U=\infty minor" );
errorbar(time_n_err, n_err, err_n, 'sqk', ...
    'MarkerFaceColor', 'k','MarkerSize',1, 'LineWidth', 1);
legend([p1 p2 p3 p4 p5 p6 p7],'Location','best');
xlabel("\it t\rm, \mus");
ylabel("\it n_{\rm O_2}\rm, mmol/m^3");
xlim([0 max(time_n_exp)]);
ylim([0 max(data_main(i_vibr,2,var,rel).nO2*1e3) + 4]);
hold off
grid minor
end

end

if P_plot
if var==1
time_p_exp=data_experiment(var).p(:,1);
p_exp=data_experiment(var).p(:,2);
end
%Pressure
if var==1
PPP(:,1)=[75 53 30 130 97 33 57 41 34];
PPP(:,2)=[0.12 0.30 0.36 0.09 0.23 0.29 0.10 0.15 0.24];
figure
hold on
plot(time_p_exp, p_exp, 'DisplayName',"p - raw data");
plot(time_p_exp(16:end), time_p_exp(16:end)*PPP(1,2)+PPP(1,1), ...
    'k-','LineWidth', 1.5, 'DisplayName',"p - interpolated data");
plot(data_main(i_vibr,2,1,rel).time, data_main(i_vibr,2,1,rel).p, ...
    'r-', 'LineWidth', 1.5, 'DisplayName', "p - U=D/6k SSH");
plot(data_main(i_vibr,3,1,rel).time, data_main(i_vibr,3,1,rel).p, ...
    'b-', 'LineWidth', 1.5, 'DisplayName', "p - U=3T SSH");
plot(data_main(i_vibr,4,1,rel).time, data_main(i_vibr,4,1,rel).p, ...
    'm-', 'LineWidth', 1.5, 'DisplayName', "p - U=\infty SSH");
legend('Location','best');
xlabel("t, \mus");
ylabel("p, Torr");
xlim([-20 100]);
ylim([-10 100]);
hold off;
grid minor;
end
end


%%
TEMP=[8110 10470 11410 7840 10310 13830 6230 7940 9560];
exp_err=[0.043 0.031 0.054 0.045 0.12 0.16 0.053 0.062 0.085];

if Nerror_plot
figure
hold on
for var=1:9

%separation of data and errors
j=1;
while (data_experiment(var).n(j,3)~=0)
j=j+1;
end
time_n_exp=data_experiment(var).n(1:j-1,1);
tlim_n=max(time_n_exp);
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

i_vibr=1; rel=2;
s1=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH D/6k
s2=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH 3T
s3=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH inf
i_vibr=2;
s4=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO D/6k
s5=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO 3T
s6=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO inf

rel=1;
i_vibr=1;
s7=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH D/6k
s8=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH 3T
s9=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %SSH inf
i_vibr=2;
s10=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO D/6k
s11=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO 3T
s12=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).nO2*1e3, ...
    min(time_n_exp):0.1:max(time_n_exp)); %FHO inf


s_exp=makima(time_n_exp, n_exp, min(time_n_exp):0.1:max(time_n_exp)); %EXP

err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
err7=mean(abs(s7-s_exp)./s_exp);
err8=mean(abs(s8-s_exp)./s_exp);
err9=mean(abs(s9-s_exp)./s_exp);
err10=mean(abs(s10-s_exp)./s_exp);
err11=mean(abs(s11-s_exp)./s_exp);
err12=mean(abs(s12-s_exp)./s_exp);

p1=plot(TEMP(var), err1, 'ok');
p2=plot(TEMP(var), err4, 'o','color', [0 0.6 0]);
p3=plot(TEMP(var), err7, 'ob');
p4=plot(TEMP(var), err10,'o','color', [0.9 0 0]);
p5=plot(TEMP(var), err2, 'sk');
p6=plot(TEMP(var), err5, 's','color',[0 0.6 0]);
p7=plot(TEMP(var), err8, 'sb');
p8=plot(TEMP(var), err11, 's','color', [0.9 0 0]);
p9=plot(TEMP(var), err3, 'dk');
p10=plot(TEMP(var), err6,'d','color',[0 0.6 0]);
p11=plot(TEMP(var), err9, 'db');
p12=plot(TEMP(var), err12, 'd','color', [0.9 0 0]);
end
set(gca, 'FontName', 'Times New Roman');
[TEMP, I]=sort(TEMP);
exp_err=exp_err(I);
p13=plot(TEMP, exp_err, 'k--');
hold off
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13],"SSH, \it U=D/6k", ...
    "FHO,\it U=D/6k", "SSH,\it U=D/6k,\rm rel off", ...
    "FHO,\it U=D/6k,\rm rel off", "SSH,\it U=3T", "FHO,\it U=3T", ...
    "SSH,\it U=3T\rm, rel off", "FHO,\it U=3T\rm, rel off", ...
    "SSH,\it U=\infty", "FHO,\it U=\infty", ...
    "SSH,\it U=\infty\rm, rel off", "FHO,\it U=\infty\rm, rel off", ...
    "mean exp. error", 'Location','eastoutside');
xlabel("\it T^0 \rm behind reflected SW, K");
xticks([6000 7000 8000 9000 10000 11000 12000 13000 14000]);
ylabel("\Delta \it n_{\rm O_2}/n_{\rm O_2}");
%title("O2/Ar mean \Delta n_{O_2}/n_{O_2}");
xlim([6000 14000]);
ylim([0 0.4]);
grid minor
end

%%
if TVerror_plot
%same as Nerror_plot
figure
set(gca,'FontSize', 12);
hold on
for var=1:9

j=1;
while (data_experiment(var).Tv(j,3)~=0)
j=j+1;
end
time_Tv_exp=data_experiment(var).Tv(1:j-1,1);
Tv_exp=data_experiment(var).Tv(1:j-1,2);
[time_Tv_exp, I]=sort(time_Tv_exp);
Tv_exp=Tv_exp(I);

k=1;
for i=1:length(Tv_exp)
    if time_Tv_exp(i)<0
        k=i;
    end
end
time_Tv_exp=time_Tv_exp(k+1:end);
Tv_exp=Tv_exp(k+1:end);

i=2;
while i~=length(time_Tv_exp)
    if time_Tv_exp(i)==time_Tv_exp(i-1)
        time_Tv_exp=[time_Tv_exp(1:i-1),
            time_Tv_exp(i+1:end)];
        Tv_exp=[Tv_exp(1:i-1),
            Tv_exp(i+1:end)];
    end
    i=i+1;
end

i_vibr=1; rel=2;
s1=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH D/6k
s2=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH 3T
s3=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH inf
i_vibr=2;
s4=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO D/6k
s5=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO 3T
s6=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO inf

rel=1;
i_vibr=1;
s7=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH D/6k
s8=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH 3T
s9=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %SSH inf
i_vibr=2;
s10=spline(data_main(i_vibr,2,var,rel).time, ...
    data_main(i_vibr,2,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO D/6k
s11=spline(data_main(i_vibr,3,var,rel).time, ...
    data_main(i_vibr,3,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO 3T
s12=spline(data_main(i_vibr,4,var,rel).time, ...
    data_main(i_vibr,4,var,rel).Tv, ...
    min(time_Tv_exp):0.1:max(time_Tv_exp)); %FHO inf


s_exp=makima(time_Tv_exp, Tv_exp, min(time_Tv_exp):0.1:max(time_Tv_exp));

err1=mean(abs(s1-s_exp)./s_exp);
err2=mean(abs(s2-s_exp)./s_exp);
err3=mean(abs(s3-s_exp)./s_exp);
err4=mean(abs(s4-s_exp)./s_exp);
err5=mean(abs(s5-s_exp)./s_exp);
err6=mean(abs(s6-s_exp)./s_exp);
err7=mean(abs(s7-s_exp)./s_exp);
err8=mean(abs(s8-s_exp)./s_exp);
err9=mean(abs(s9-s_exp)./s_exp);
err10=mean(abs(s10-s_exp)./s_exp);
err11=mean(abs(s11-s_exp)./s_exp);
err12=mean(abs(s12-s_exp)./s_exp);

p1=plot(TEMP(var), err1, 'ok');
p2=plot(TEMP(var), err4, 'o','color', [0 0.6 0]);
p3=plot(TEMP(var), err7, 'ob');
p4=plot(TEMP(var), err10,'o','color', [0.9 0 0]);
p5=plot(TEMP(var), err2, 'sk');
p6=plot(TEMP(var), err5, 's','color',[0 0.6 0]);
p7=plot(TEMP(var), err8, 'sb');
p8=plot(TEMP(var), err11, 's','color', [0.9 0 0]);
p9=plot(TEMP(var), err3, 'dk');
p10=plot(TEMP(var), err6,'d','color',[0 0.6 0]);
p11=plot(TEMP(var), err9, 'db');
p12=plot(TEMP(var), err12, 'd','color', [0.9 0 0]);
end
set(gca, 'FontName', 'Times New Roman');
[TEMP, I]=sort(TEMP);
exp_err=exp_err(I);
p13=plot(TEMP, exp_err, 'k--');
hold off
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13], "SSH, D/6k", ...
    "FHO, D/6k", "SSH, D/6k, rel off", "FHO, D/6k, rel off"...
 , "SSH, 3T", "FHO, 3T", "SSH, 3T, rel off", "FHO, 3T, rel off"...
 , "SSH, inf", "FHO, inf", "SSH, inf, rel off", "FHO, inf, rel off", ...
 "mean exp. error", 'Location','eastoutside');
xlabel("T^0 behind reflected SW, K");
xticks([6000 7000 8000 9000 10000 11000 12000 13000 14000]);
ylabel("\DeltaT_v/T_v");
%title("O2/Ar mean \DeltaT_v/T_v");
xlim([6000 14000]);
ylim([0 0.6]);
grid minor
end