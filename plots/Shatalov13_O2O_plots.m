%% переделываем графики по Шаталову для STS-Simple-Solver
% 30.08.2023 Maksim Melnik

%% попытка найти лучшую модель диссоциации для кислорода для правок по 
    % статье АА, май 2020
% содержание:
% 1) How to use?
% 2) loading exp. data
% 3) считаем максимальное относительное отклонение для Tv в %
% 4) выводим данные по максимальным максимальным отклонениям
% 5) рисуем сравнение Tv для разных моделей
% 12) рисуем график отклонения от эксперимента
%% How to use?
% 1) Choose desired test cases in examples/O2O_SW_Shatalov.m
%       To simulate the whole set of combinations for models use i_U = 2:4
%       i_vibr = 1:3. The whole set of test cases is i_ini = 1:5.
% 2) Run examples/O2O_SW_Shatalov.m. e.g. res_Shatalov = O2O_SW_Shatalov;
% 3) Choose below out = *your output*, e.g. out = res_Shatalov;
% 4) Choose models in the block 5 for plots and run this block. You can
%       turn off the benchmark.
%% 2) loading
res_Shatalov = O2O_SW_Shatalov;
load("../data/for comparison/shatalov.mat")
out = res_Shatalov;
str_vibr = ["SSH", "FHO", "FHO-FR-reg"];
str_U    = ["Savelev", "\itD/\rm6\itk", "3\itT", "\infty"];
pl_vibr  = [":",   "--",  "-"];
pl_U     = [0.5 0.5 0
            0   0.6 0
            0.9 0   0
            0   0   0.6
            0   0.5 0.5];
%% 3 считаем максимальное относительное отклонение для Tv в %
dTv=zeros(5,5,5);
for model_vibr = 2:3 % 1 - SSH, 2 - FHO, 3 - FHO-FR
for model_U = 2:4
  for ind_ini=1:5
      interP=shatalov(ind_ini).Tv(:,1);
      temp_dT=spline(out(ind_ini, model_vibr, model_U).time_mus, ...
                            out(ind_ini, model_vibr, model_U).Tv, interP);
      del=(temp_dT-shatalov(ind_ini).Tv(:,2))./shatalov(ind_ini).Tv(:,2);
      if ind_ini==2
          del(2)=0;
      end
      temp_dT=max(abs(del));
      % temp_dT = mean(abs(del));
      dTv(ind_ini, model_vibr, model_U)=temp_dT;
      clear temp_dT;
  end
end
end
dTv=dTv*100;
%% 4 выводим данные по максимальным максимальным отклонениям
disp('--------------------')
for model_vibr = 1:3
	for model_U = 2:4
		disp(['model_vibr=' num2str(model_vibr) ', model_U=' num2str(model_U)...
            ', max dT=' num2str(max(dTv(:, model_vibr, model_U)))])
	end
end
%% 5 рисуем сравнение Tv для разных моделей
FontSize=16;
FigureSize=[0 0 700 600];
FigureSize=[20 50 600 560];
for ind_ini = 2 % [2 4]   % chossing experimental test cases
 figure('Units', 'pixels', 'OuterPosition', FigureSize);
 set(gca, 'FontName', 'Times New Roman');
 hold on;%grid on;box on;
 grid minor
  % plotting an experimental data
 kkk=zeros(length(shatalov(ind_ini).Tv(:,2)),1)+0.1;
 if ind_ini<4
    kkk=kkk+0.05;
    for jj=1:length(shatalov(ind_ini).Tv(:,2))
       if shatalov(ind_ini).Tv(jj,1)<=0.2
          kkk(jj)=0.25;
       end
    end
 end
 [~, imax]=max(shatalov(ind_ini).Tv(:,2));
 kkk2=kkk;
 mltpl=0.2;
 kkk2(1:imax)=kkk2(1:imax)+abs((1:imax)-imax)'/(imax-1)*mltpl+0.1;
 sec_i=round((length(kkk)-imax)/4*3);
 kkk2(imax+1:imax+sec_i)=kkk2(imax+1:imax+sec_i)+...
    abs((1:sec_i)-sec_i)'/(sec_i-1)*0.1;
 kkk3=kkk*0;
 if ind_ini > 3
    kkk3(end-2)=kkk(end-2)+0.02;
    kkk3(end-1)=kkk(end-1)+0.03;
    kkk3(end)=kkk(end)+0.06;
 end
 p_e1=errorbar(shatalov(ind_ini).Tv(:,1),shatalov(ind_ini).Tv(:,2),...
    kkk2.*shatalov(ind_ini).Tv(:,2), 'ksq', 'MarkerFaceColor', 'k',...
     'CapSize',4.5, 'LineWidth',0.6, ...
     'DisplayName', "ISW, exp");
  % plotting the benchmark
 model_vibr_bm = 2; % benchmark vibrational model
 model_U_bm = 3;    % benchmark U parameter
 % p_bm = plot(out(ind_ini, model_vibr_bm, model_U_bm).time_mus, ...
 %                        out(ind_ini, model_vibr_bm, model_U_bm).Tv,...
 %         pl_vibr(model_vibr_bm), 'Color', pl_U(model_U_bm, :),'LineWidth',1.5, 'DisplayName', ...
 %         "Benchmark: " + str_vibr(model_vibr_bm) + ", U = " + str_U(model_U_bm));

 for model_vibr=2:3 %1:3  % choosing vibrational models to plot
  for model_U = 2:4     % choosing U parameter for plotting
      % the plot
    p3=plot(out(ind_ini, model_vibr, model_U).time_mus, out(ind_ini, model_vibr, model_U).Tv,...
         pl_vibr(model_vibr), 'Color', pl_U(model_U, :),'LineWidth',2, 'DisplayName', ...
         str_vibr(model_vibr) + ", " + str_U(model_U));
      % some plot's parameters
    set(gca, 'FontSize', FontSize);
    ylabel('\it T_{\rmv}^{\rm O_2}, \rmK')%,'FontSize',FontSize);
    xlabel('\itt, \rmмкс')%,'FontSize',FontSize);
    xlim([min(shatalov(ind_ini).Tv(:,1))-0.1 max(shatalov(ind_ini).Tv(:,1))+0.2]);
    set(gca, 'TickDir', 'In',... чёрточки внутри, out -- снаружи
        'LineWidth', 1,... толщина окантовки
        'GridAlpha', 0.2) % прозрачность сетки
  end
 end
 %title("T = " + fix(out(ind_ini, model_vibr, model_U).T(1)) + " K")
 legend('Location', 'best');
end
%% рисуем график отклонения от эксперимента
numM = [13.5, 12.3, 12, 10.2, 9.7];
% FigureSize=[0 0 600 560];
clr=[   0 0 0
        1 0 0
        0 0 1
        0 0.5 0
        1 0 1];
clr = pl_U;
mark=["sq" "o" "^"];
figure('Units', 'pixels', 'OuterPosition', FigureSize);
 set(gca, 'FontName', 'Times New Roman');
hold on;grid on;%box on;
 % grid minor
markersize=9;
markerwidth=1.5;
modV_benchmark = 2;
modU_benchmark = 4;
plot(numM, reshape(dTv(:, modV_benchmark, modU_benchmark), [], 1), ...
    mark(modV_benchmark), ...
        'Color', clr(modU_benchmark, :), 'MarkerSize', markersize, ...
        'LineWidth', markerwidth);
for modV = 3
    plot(numM(1:end), reshape(dTv(:, modV, 2), [], 1), mark(modV), ...
        'Color', clr(2, :), 'MarkerSize', markersize, ...
        'LineWidth', markerwidth);
    plot(numM(1:end), reshape(dTv(:, modV, 3), [], 1), mark(modV),...
        'Color', clr(3,:), 'MarkerSize', markersize, ...
        'LineWidth', markerwidth);
    plot(numM(1:end), reshape(dTv(:, modV, 4), [], 1), mark(modV),...
        'Color', clr(4,:), 'MarkerSize', markersize, ...
        'LineWidth', markerwidth);
%     plot(numM(1:end),reshape(dTv(modV,1,:),[],1),mark(modV),...
%         'Color',clr(1,:), 'MarkerSize', markersize, ...
%         'LineWidth', markerwidth);
end
legend('FHO, \itU\rm = \infty','FHO-FR-reg, \itU = D/\rm6\itk', ...
    'FHO-FR-reg, \itU = \rm3\itT', ...
    'FHO-FR-reg, \itU\rm = \infty',...
    ...'SSH-Savelev', 'FHO-Pogosbekyan','FHO-\infty','FHO-3T','FHO-D/6k',...
    ...'FHO-Savelev', 'FHO2-Pogosbekyan','FHO2-\infty','FHO2-3T',...
    ...'FHO2-D/6k','FHO2-Savelev',...
                             'Location','best','FontSize',FontSize);
set(gca, 'LineWidth', 1,... толщина окантовки
         'GridAlpha', 0.2) % прозрачность сетки
xlim([min(numM(1:end))-0.2 max(numM)+0.2]);
ylim([0 45])
set(gca, 'FontSize', FontSize);
ylabel('макс. отклонение    \itT\rm_v, %')
xlabel('\itM\rm_0')
yticks([0 10 15 25 40])
xticks([numM(5) numM(4) 11 numM(3) numM(2) numM(1)])
% 
% % добавочка для графика средних
% ylim([0 30])
% ylabel('mean \delta_{O_2} (%)')