% tc=tc4;
% tc=res_A;
% tc=res_MT;
% tcwd=res_MTwd; % without dissociation
% tc1 -- новые статвеса для дет баланса (пофиксить экспоненты)
% tc2 -- старые статвеса для дет баланса
% tc3 -- старые статвеса + исправленная экспонента с энергиями электр.сост
% tc5 -- без электронных состояний CO
%
% How to use?
% - load res_all
% - load num_VARIABLES from begining of main code (i_dis>=2)
% - choose params
% - run code
% To draw plots for paper and report for CMMASS comparing diss models 
%   use 14.
% To draw plots for paper and report for CMMASS comparing U use 15.
% Overall:
% 1. Temperature plot
% 2. rho CO(X) по старому дет балансу
% 3. rhoCOa/rho
% 4. rhoCOA/rho
% 5. n_C/n
% 6. fig 2 (Aliat 2003), vibr levels
% 7. fig 7 (Aliat 2003), with or without diss for COX, i=40
% 8. fig 10 (Aliat 2003) for mass fraction COa with and without diss
% 9. fig 15 (Aliat 2003) to check COA, COa, C and O for full scheme
% 10. seaching for the induction time
% 11. draw induction time
% 12. draw nC for different diss models
% 13. draw nC just for one case
% 14. draw nC for different diss models for format res(i_ini, i_dis, i_U)
% 15. draw nC for different U for format res(i_ini, i_dis, i_U)
% 16. draw nCO just for one case
% 17. draw tau Appleton
% 18. draw nCO(1-alpha) like Fairbairn fig 10
% 18.1. draw nCO like Fairbairn fig 9
% 18.2. draw nCO(1-alpha) like Fairbairn fig 11
% 19. draw nCO like Appleton fig 3
% 20. draw nC like Appleton fig 3 for 3 diss models
% 22. five test cases from Fairbairn's paper
% 24. universal drawing results script
% 25. draw nC2
% 26. draw nC2 for Fairbairn 8f plot 11
% 27. draw nC2 for Fairbairn 8e plot 12
load res_all
%% 1. temperature
plot(tc.temp(:,1), tc.temp(:,1+num_T), 'linewidth', 1.5)
xlim([0 0.5])
ylim([7000 12000])
%% 2. rho CO(X) по старому дет балансу
x0p001=find(tc.temp(:,1)<0.001, 1, 'last');
% rho incorrect
% rho=sum(tc.temp(x0p001,2:CO.num_vibr_levels(3)+1))*CO.mass...
%     +tc.temp(x0p001, num_C+1)*C.mass+tc.temp(x0p001, num_O+1)*O.mass...
%    +tc.temp(x0p001, num_C2+1)*C2.mass+tc.temp(x0p001, num_Ar+1)*Ar.mass;
semilogy(0:CO.num_vibr_levels(1)-1, ...
    tc.temp(x0p001,2:CO.num_vibr_levels(1)+1)*CO.mass/rho,...
                                                        'linewidth', 1.5)
xlim([0 40])
ylim([1e-6 1])
%% 3. rhoCOa/rho
% rhox probably incorrect
% rhox=sum(tc.temp(:,2:1+num_COa+CO.num_vibr_levels(2)-1), 2)*CO.mass...
%     +tc.temp(:, num_C+1)*C.mass+tc.temp(:, num_O+1)*O.mass...
%     +tc.temp(:, num_C2+1)*C2.mass+tc.temp(:, num_Ar+1)*Ar.mass;
semilogy(tc.temp(:, 1), ...
 sum(tc.temp(:,1+num_COa:1+num_COa+CO.num_vibr_levels(2)-1), 2)*CO.mass...
                                                ./rhox,'linewidth', 1.5)
xlim([0 0.3])
ylim([1e-5 1e-3])
%% 4. rhoCOA/rho
rhox=sum(tc.temp(:,2:1+num_COA+CO.num_vibr_levels(3)-1), 2)*CO.mass...
    +tc.temp(:, num_C+1)*C.mass+tc.temp(:, num_O+1)*O.mass...
    +tc.temp(:, num_C2+1)*C2.mass+tc.temp(:, num_Ar+1)*Ar.mass;
semilogy(tc.temp(:, 1), ...
 sum(tc.temp(:,1+num_COA:1+num_COA+CO.num_vibr_levels(3)-1), 2)*CO.mass...
                                                ./rhox,'linewidth', 1.5)
% xlim([0 0.3])
% ylim([1e-8 1e-4])
% ylim([1e-8 1e-6])
ylim([1e-12 1e-7])
% title('COA without diss from COA and COa')
%% 5. n_C/n
v0=5200;
time_ms=tc.temp(:, 1)/v0*1e6;
plot(time_ms, ...
  tc.temp(:,1+num_C)./sum(tc.temp(:,1:1+num_Ar), 2),'r','linewidth', 1.5)
% semilogx(time_ms, ...
%      tc.temp(:,1+num_C)./sum(tc.temp(:,1:1+num_Ar), 2),'linewidth', 1.5)
xlim([1e-2 1e3])
xlim([1e1 1e4])
xlim([1e1 1e3])
ylim([0 1e-1])
ylim([0 6e-3])

%% draw plots for comparison with Aliat 2003
%% 6. fig 2, vibr levels
ind_i=[0 1 5 10 15 20 30 40];
rho=sum(tc.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tc.temp(:, 1+num_C)*C.mass+tc.temp(:, 1+num_O)*O.mass...
   +tc.temp(:, 1+num_C2)*C2.mass+tc.temp(:, 1+num_Ar)*Ar.mass;
loglog(tc.temp(:,1)*1e2, tc.temp(:, 2)*CO.mass./rho)
hold on
for i=ind_i(2:end)
 loglog(tc.temp(:,1)*1e2, tc.temp(:, 2+i)*CO.mass./rho)
end
ylim([1e-8 1e0])
xlim([1e-3 1e1])
%% 7. fig 7, with or without diss for COX, i=40

rho=sum(tc.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tc.temp(:, 1+num_C)*C.mass+tc.temp(:, 1+num_O)*O.mass...
   +tc.temp(:, 1+num_C2)*C2.mass+tc.temp(:, 1+num_Ar)*Ar.mass;
rhowd=sum(tcwd.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tcwd.temp(:, 1+num_C)*C.mass+tcwd.temp(:, 1+num_O)*O.mass...
   +tcwd.temp(:, 1+num_C2)*C2.mass+tcwd.temp(:, 1+num_Ar)*Ar.mass;
i=40;
semilogx(tc.temp(:,1)*1e2, tc.temp(:,2+i)*CO.mass./rho, ...
    tcwd.temp(:,1)*1e2, tcwd.temp(:,2+i)*CO.mass./rhowd)
legend('VT+VE+diss', 'VT+VE', 'location', 'best')
xlim([1e-1 1e2])
ylim([0 3e-5])
%% 8. fig 10 for mass fraction COa with and without diss
rho=sum(tc.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tc.temp(:, 1+num_C)*C.mass+tc.temp(:, 1+num_O)*O.mass...
   +tc.temp(:, 1+num_C2)*C2.mass+tc.temp(:, 1+num_Ar)*Ar.mass;
rhowd=sum(tcwd.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tcwd.temp(:, 1+num_C)*C.mass+tcwd.temp(:, 1+num_O)*O.mass...
   +tcwd.temp(:, 1+num_C2)*C2.mass+tcwd.temp(:, 1+num_Ar)*Ar.mass;
COa_f=sum(tc.temp(:, 1+num_COa:num_COA), 2)*CO.mass./rho;
COa_fwd=sum(tcwd.temp(:, 1+num_COa:num_COA), 2)*CO.mass./rhowd;
semilogy(tc.temp(:,1)*1e2, COa_f, ...
    tcwd.temp(:,1)*1e2, COa_fwd)
legend('VT+VE+diss', 'VT+VE', 'location', 'best')
ylim([1e-4 1e-2])
xlim([0 10])
%% 9. fig 15 to check COA, COa, C and O for full scheme
rho=sum(tc.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
   +tc.temp(:, 1+num_C)*C.mass+tc.temp(:, 1+num_O)*O.mass...
   +tc.temp(:, 1+num_C2)*C2.mass+tc.temp(:, 1+num_Ar)*Ar.mass;
COa_f=sum(tc.temp(:, 1+num_COa:num_COA), 2)*CO.mass./rho;
COA_f=sum(tc.temp(:, 1+num_COA:1+num_COA+CO.num_vibr_levels(3)-1), 2)...
    *CO.mass./rho;
C_f=tc.temp(:, 1+num_C)*C.mass./rho;
O_f=tc.temp(:, 1+num_O)*O.mass./rho;
semilogy(tc.temp(:,1)*1e2, O_f, tc.temp(:,1)*1e2, C_f, ...
    tc.temp(:,1)*1e2, COa_f, tc.temp(:,1)*1e2, COA_f)
legend('O', 'C', 'COa', 'COA', 'location', 'best')
ylim([1e-10 1e0])
xlim([0 30])

%% 10. here seaching for the induction time
v0=5200; disp('Aliats testcase')
% v0=2860;
time_ms=tc.temp(:, 1)/v0*1e6;
plot(time_ms, tc.temp(:,1+num_C),'r','linewidth', 1.5)
xlim([0 2e1])
% ylim([0 5e-6])
xticks(1:0.5:20)
xtickangle(90)
%% 11. draw induction time
tau=1.9;

%% 12. draw nC for different diss models
% v0=5200; disp('Aliats testcase')
v0=2860;
% v0=2735;
time_ms_MT=res_MT.temp(:, 1)/v0*1e6;
time_ms_A=res_A.temp(:, 1)/v0*1e6;
time_ms_woE=res_woE.temp(:, 1)/v0*1e6;
figure
plot(time_ms_MT, res_MT.temp(:, 1+num_C), ...
    time_ms_A, res_A.temp(:,1+num_C),...
    time_ms_woE, res_woE.temp(:, 1+num_C),...
    'linewidth', 1.5)
xlim([1e-2 1e3]) % Aliat
% xlim([1e1 1e5])
xlim([0 1e1])
legend('COX, COa, COA MT', 'COX, COa, COA Aliat', 'COX MT', ...
    'location', 'best')
% legend('D/6k', '3T', 'Inf', 'location', 'best')
%% 13. draw nC just for one case
% res_all(7,1,2)=res_7(7,1,2);
i_ini=5;
i_dis=4;
i_U=3;
plot(res_all(i_ini, i_dis, i_U).temp(:, end),...
            res_all(i_ini, i_dis, i_U).temp(:, 1+num_C), 'linewidth', 1.5)
% ylim([0 6e13])
xlim([0 35])
%% 14. draw nC for different diss models for format res(i_ini, i_dis, i_U)
%, ...
%                             'color', [126 47 142]/255, 'linewidth', 1.5);
i_ini=2;
% i_dis=3;
i_U=3;
i_plot=1;
figure
axesH=axes;gbuild_fun(res_all(i_ini, 1, i_U).temp(:, end), ...
    res_all(i_ini, 1, i_U).temp(:, 1+num_C), [33 115 183]/255, i_plot);
hold on
gbuild_fun(res_all(i_ini, 2, i_U).temp(:, end), ...
    res_all(i_ini, 2, i_U).temp(:, 1+num_C), [217 83 25]/255, i_plot);
gbuild_fun(res_all(i_ini, 3, i_U).temp(:, end), ...
    res_all(i_ini, 3, i_U).temp(:, 1+num_C), [255 170 0]/255, i_plot);
gbuild_fun(res_all(i_ini, 4, i_U).temp(:, end), ...
    res_all(i_ini, 4, i_U).temp(:, 1+num_C), [126 47 142]/255, i_plot);
xlim([7e-2 3e2])    % Aliat(3), logarithmic
% ylim([0 1.7e23])
% xlim([1e1 1e8])
% xlim([1e3 5e7])     % Appleton, 5, logarithmic
% ylim([0 4.5e19])
% xlim([0 50])
xlim([0 35])        % Aliat(3), Appleton(5), linear
% xlim([0 0.6])        % Fairbairn(2)
FontSize=16;
set(gca, 'LineWidth', 1,...         толщина окантовки
         'GridAlpha', 0.2,...
         'FontSize', FontSize,...   % прозрачность сетки
         'TickDir', 'Out',...       % чёрточки наружу
         'Box', 'off',...           % сверху-справа нет шкалы
         'FontName', 'Times'...    % шрифт
         ,'XMinorTick', 'on'...
         ...,'XTick', [1e-1 1e0 1e1 1e2] ... Aliat
          ...,'XTick', [1e3 1e4 1e5 1e6 1e7] ... Appleton
         ,'XTick', 0:10:30 ...
    )
axesH.XAxis.MinorTickValues = 5:10:35;
legend('No electronic excitation',...Kinetic scheme 1',...
    ...'CO($X^1\Sigma^+$) M-T', ...
    'Kinetic scheme 2',...'CO($X^1\Sigma^+, a^3\Pi, A^1\Pi$) M-T', ...
    'Kinetic scheme 3',...'CO($X^1\Sigma^+, a^3\Pi, A^1\Pi$) Aliat',...
    'Kinetic scheme 4',...CO(a), CO(A) Aliat', ...
    'location', 'best'..., 'Interpreter','latex'...
    )
legend('boxoff')
xlabel('t, мкс')                                    % in RU
ylabel('n_{C}, м^{-3}')
xlabel('$t, \mu$s', 'Interpreter','latex')          % in EN
ylabel('$n_{C}$, m$^{-3}$', 'Interpreter','latex')
% legend('D/6k', '3T', 'Inf', 'location', 'best')
%% 15. draw nC for different U for format res(i_ini, i_dis, i_U)
i_ini=2;
i_dis=3;
% i_U=2;
figure
axesH=axes;
plot(res_all(i_ini, i_dis, 2).temp(:, end),...
                    res_all(i_ini, i_dis, 2).temp(:, 1+num_C),...
    res_all(i_ini, i_dis, 3).temp(:, end),...
                    res_all(i_ini, i_dis, 3).temp(:, 1+num_C),...
    res_all(i_ini, i_dis, 4).temp(:,end),...
            res_all(i_ini, i_dis, 4).temp(:, 1+num_C), 'linewidth', 1.5)
xlim([5e-2 3e2]) % Aliat
xlim([5e1 1e5])
% xlim([0 35])
xlim([0 20])    % Fairbairn(2), linear
FontSize=16;
set(gca, 'LineWidth', 1,...           толщина окантовки
         ...'GridAlpha', 0.2,...
         'FontSize', FontSize,...   % прозрачность сетки
         'TickDir', 'Out',...       % чёрточки наружу
         'Box', 'off',...           % сверху-справа нет шкалы
         'FontName', 'Times'...       % шрифт
         ,'XMinorTick', 'on'...
         ,'XTick', 0:10:30 ...
         )
axesH.XAxis.MinorTickValues = 5:10:35;
xlabel('t, мкс')                               % in RU
ylabel('n_{C}, м^{-3}')
xlabel('$t, \mu$s', 'Interpreter','latex')          % in EN
ylabel('$n_{C}$, m$^{-3}$', 'Interpreter','latex')
legend('$U=D/6k$', '$U=3T$', '$U=\infty$', 'location', 'best', ...
    'Interpreter','latex')
legend('boxoff')
%% 16. draw nCO just for one case
res=res_A;
v0=2735;
time_ms=res.temp(:, 1)/v0*1e6;
nCO=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
plot(time_ms, nCO/1e6, ...
    'linewidth', 1.5)
% ylim([0 6e13])
xlim([0 800])
%% 17. draw tau Appleton
res=res_AAp3;
T=res.temp(end, 1+num_T);
nCO=sum(res.temp(1, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
nM=nCO+res.temp(1, 1+num_Ar);
% p=10*Torr;
% n=p/V_K/T;
% nCO=n*0.01;
% nM=n;
TAp=500:500:17000;
fAp=4.97e-19*exp(-132400./TAp);
tau2CAr=1./fAp;
tau2=tau2CAr/(nCO*1e-6)/(nM*1e-6);
tau_2=sqrt(tau2)*1e6;
semilogy(1./TAp*1e4, fAp.^0.5)
xlim([0.6 1.8])
tau=sqrt(1/(4.97e-19*exp(-132400/T)/1e12*nM*nCO))*1e6;
% hold on
% plot([tau tau], [1e13 1.2e13])
%% 18. draw nCO(1-alpha) like Fairbairn fig 10
% !rescale like in Appleton paper
res=res_all(2,4,3);
v0=2860; disp('Fairbairn 8a')
nCO=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
nC=res.temp(:, 1+num_C);
nCOalpha=nCO./(nCO+nC);
time_ms=res.temp(:, 1)/v0*1e6/4;
plot(time_ms, nCOalpha)
xlim([0 40])
ylim([0.5 1])
%% 18.1. draw nCO like Fairbairn fig 9
% !rescale like in Appleton paper
res=res_all(2,4,3);
v0=2860; disp('Fairbairn 8a')
nCO=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
nC=res.temp(:, 1+num_C);
% nCOalpha=nCO./(nCO+nC);
time_ms=res.temp(:, 1)/v0*1e6/4;  % musec, but like in fig 9
plot(time_ms, nCO)
xlim([0 40])
% ylim([0.5 1])
%% 18.2 draw nCO(1-alpha) like Fairbairn fig 11
% !rescale like in Appleton paper
res=res_all(8,1,3);
v0=3260; disp('Fairbairn 8f')
nCO=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
nC=res.temp(:, 1+num_C);
nCOalpha=nCO./(nCO+nC);
time_ms=res.temp(:, 1)/v0*1e6/4;
plot(time_ms, nCOalpha)
xlim([0 20])
% ylim([0.5 1])
%% 19. draw nCO like Appleton fig 3
res=res_AAp3;
% res=res_MTAp3;
% res=res_woEAp3;
v0=2090;    disp('Appleton fig3')
nCO=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)-1), 2);
time_ms=res.temp(:, 1)/v0*1e6;
% rho=sum(res.temp(:, 2:1+num_COA+CO.num_vibr_levels(3)), 2)*CO.mass...
%    +res.temp(:, 1+num_C)*C.mass+res.temp(:, 1+num_O)*O.mass...
%    +res.temp(:, 1+num_C2)*C2.mass+res.temp(:, 1+num_Ar)*Ar.mass;
semilogy(time_ms, nCO)
xlim([0 35])
% ylim([max(nCO)/4 max(nCO)/4*10])
figure
nC=res.temp(:, 1+num_C);
time_ms=res.temp(:,end);
plot(time_ms, nC)
xlim([0 35])
%% 20. draw nC like Appleton fig 3 for 3 diss models
% plot(res_AAp3.temp(:,end), res_AAp3.temp(:,1+num_C),...
%     res_MTAp3.temp(:,end), res_Ap3.temp(:,end)
%% 21. 
res_AAp3;    % Appleton fig 3, diss Aliat
res_A8a;     % Fairbairn fig 8a, diss Aliat
%% 22. Я взял пять случаев из статьи Фейрбеирна
%             %   f       p, Torr     v0, m/sec   T, K  
%     init_c=[    0.01    10          2860        7600    % 8a
%                 0.01    20          2820        7400    % 8b
%                 0.02    3           3000        8300    % 8c
%                 0.001   10          3080        8800    % 8e
%                 0.1     5.4         3260        9000    % 8f ];
%% 23. draw 0th and 1st states to check incorrect Tv values in my 
...calculations
temp=res(7, 3, 2).temp;
time=temp(:,end);
time=temp(:,1);
nCO0=temp(:,1+1);
nCO1=temp(:,1+2);
Tv=temp(:,end-1);
figure
semilogx(time, Tv)
% loglog(time, nCO0)
% loglog(time, nCO1)
%% 24. draw nC for different diss models for format res(i_ini, i_dis, i_U)
str_dis=["MT, ", "MT with exc, ", "Aliat, ", "Savelev, "];
str_U=["Savelev ", "D/6k", "3T", "inf"];
for i_ini=2:2
    figure
legend_str='';
 for i_dis=1:1
  for i_U=3
semilogx(res_all(i_ini, i_dis, i_U).temp(:, end),...
    res_all(i_ini, i_dis, i_U).temp(:, 1+num_C), 'linewidth', 1.5)
hold on
name=str_dis(i_dis)+"U="+str_U(i_U);
legend_str=[legend_str name];
  end
 end
title("TC"+num2str(i_ini)+" (in code notation)")
FontSize=16;
set(gca, 'LineWidth', 1,...         толщина окантовки
         'GridAlpha', 0.2,...
         'FontSize', FontSize,...   % прозрачность сетки
         'TickDir', 'Out',...       % чёрточки наружу
         'Box', 'off',...           % сверху-справа нет шкалы
         'FontName', 'Times'...    % шрифт
         ,'XMinorTick', 'on'...
         ...,'XTick', [1e-1 1e0 1e1 1e2] ... Aliat
         ...,'XTick', [1e3 1e4 1e5 1e6 1e7] ... Appleton
         ...,'XTick', 0:10:30 ...
    )
legend(legend_str(2:end), 'location', 'best')
legend('boxoff')
% legend('No electronic excitation',...
% ...    'Kinetic scheme 2',...'CO($X^1\Sigma^+, a^3\Pi, A^1\Pi$) M-T', ...
%     'Kinetic scheme 3',...
%     'location', 'best'..., 'Interpreter','latex'...
%     )
% legend('boxoff')
xlabel('$t, \mu$s', 'Interpreter','latex')          % in EN
ylabel('$n_{C}$, m$^{-3}$', 'Interpreter','latex')
xl_arr=[5e-1 8e3;
        2e1 8e4;
        1e-2 1e3;   % Aliat
        1e2 5e5;
        1e3 5e7;     % Appleton
        3e0 2e5
        1e0 1e4
        ];
xlim(xl_arr(i_ini,:))
end
%% 25. draw nC2 for Fairbairn 8a plot
figure
res=res_all(2,1,3);
v0=2860; disp('Fairbairn 8a')
nC2=res.temp(:, 1+num_C2)/1e6*1e-13;
time_ms=res.temp(:, 1)/v0*1e6/4;    % divided by 4 because of 
plot(time_ms, nC2)                  %       Appleton remark
xlim([0 40])
% ylim([0 4.1e-13])
ylim([0 4])
title('C_2')
%% 26. draw nC2 for Fairbairn 8f plot 11
figure('Units', 'pixels', 'OuterPosition', [600 0 400 560])
res=res_all(8,1,3);
v0=3260; disp('Fairbairn 8f')
nC2=res.temp(:, 1+num_C2)/1e6*1e-13;
time_ms=res.temp(:, 1)/v0*1e6/4;    % divided by 4 because of 
plot(time_ms, nC2)                  %       Appleton remark
xlim([0 25])
% ylim([0 4.1e-13])
ylim([0 6])
title('C_2')
%% 27. draw nC2 for Fairbairn 8e plot 12
figure('Units', 'pixels', 'OuterPosition', [600 0 400 560])
res=res_all(9,1,3);
v0=3080; disp('Fairbairn 8e')
nC2=res.temp(:, 1+num_C2)/1e6*1e-13;
time_ms=res.temp(:, 1)/v0*1e6/4;    % divided by 4 because of 
plot(time_ms, nC2)                  %       Appleton remark
xlim([0 45])
% ylim([0 4.1e-13])
ylim([0 1])
title('C_2')
%% plotting functions
function gbuild_fun(x, y, c, p) 
    switch p
        case 1
            plot(x, y, 'color', c, 'linewidth', 1.5)
        case 2
            semilogx(x, y, 'color', c, 'linewidth', 1.5)
        case 3
            semilogy(x, y, 'color', c, 'linewidth', 1.5)
        case 4
            loglog(x, y, 'color', c, 'linewidth', 1.5)
    end
end