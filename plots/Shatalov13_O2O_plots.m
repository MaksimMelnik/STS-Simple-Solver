%% Plots for Shatalov (Ibraguimova) experiment [1] 
% taking into account error correction [2]. Also the recommended models
% combination [3] is included.
% [1] Ibraguimova, et al J. Chem. Phys. 139, 034317 (2013).
% [2] I.J. Wysong, et al AIP Conf. Proc. 2132, 180008 (2019).
% [3] L. Campoli, et al. Acta Astronautica 175 (2020) 493–509.
% 30.08.2023 Maksim Melnik

% contents:
% 1) How to use?
% 2) loading exp and simulated data
% 3) counting maximum deviation of Tv from exp data in %
% 4) displaying max deviation data
% 5) drawing Tv plots vs exp data
%% 1) How to use?
% 1) Choose desired test cases in examples/O2O_SW_Shatalov.m
%       To simulate the whole set of combinations for models use i_U = 2:4
%       i_vibr = 1:3. The whole set of test cases is i_ini = 1:5.
% 2) Run examples/O2O_SW_Shatalov.m. e.g. res_Shatalov = O2O_SW_Shatalov;
% 3) Choose below out = *your output*, e.g. out = res_Shatalov;
% 4) Choose models in the block 5 for plots and run this block. You can
%       turn off the benchmark.
%% 2) loading exp and simulated data
load("../data/for comparison/Ibraguimova2013.mat")     % loading exp data
    % just in case you have exp data by authors
% load("../data/for comparison/shatalov.mat")     
out      = res_Shatalov;                        % loading simulated data
str_vibr = ["SSH", "FHO", "FHO-FR"];            % vibr models
str_U    = ["Savelev", "D/6k", "3T", "\infty"]; % U values
pl_vibr  = [":",   "-",  "-."];                 % lines for plots
pl_U     = [0.5 0.5 0                           % colors of line plots
            0   0.6 0
            0.9 0   0
            0   0   0.6
            0   0.5 0.5];
%% 3) counting maximum deviation of Tv from exp data in %
dTv = zeros(5, 3, 5);
model_vibr = 2;     % 1 - SSH, 2 - FHO
for model_U = 2:4
  for ind_ini = 1:5
      interP = shatalov(ind_ini).Tv(:, 1);
      temp_dT = spline(out(ind_ini, model_vibr, model_U).time_mus, ...
                            out(ind_ini, model_vibr, model_U).Tv, interP);
      del = (temp_dT - shatalov(ind_ini).Tv(:, 2)) ./ ...
                                            shatalov(ind_ini).Tv(:, 2);
      if ind_ini == 2
          del(2) = 0;
      end
      temp_dT = max(abs(del));
      dTv(ind_ini, model_vibr, model_U) = temp_dT;
      clear temp_dT;
  end
end
dTv = dTv * 100;
%% 4) displaying max deviation data
disp('--------------------')
for model_vibr = 1:2
 for model_U = 2:4
  disp(['model_vibr = ' num2str(model_vibr) ...
        ', model_U = ' num2str(model_U) ...
                ', max dT=' num2str(max(dTv(:, model_vibr, model_U)))])
 end
end
%% 5) drawing Tv plots vs exp data
FontSize = 20;
FigureSize = [0 0 700 600];
for ind_ini = 1 %1:5 % [2 4]   % chossing experimental test cases
 figure('Units', 'pixels', 'OuterPosition', FigureSize); % window options
 hold on;           % several plots on the same plot
 grid on; box on;   % grid and frame
    % plotting an experimental data
 Tv_exp = zeros(length(shatalov(ind_ini).Tv(:, 2)), 1) + 0.1;
 if ind_ini < 4     % experimental error [1]
    Tv_exp = Tv_exp + 0.05;
    for jj = 1:length(shatalov(ind_ini).Tv(:, 2))
       if shatalov(ind_ini).Tv(jj, 1) <= 0.2
          Tv_exp(jj) = 0.25;
       end
    end
 end
 [~, imax]  = max(shatalov(ind_ini).Tv(:, 2));
 Tv_exp_2   = Tv_exp;
 mltpl      = 0.2;
 Tv_exp_2(1:imax) = Tv_exp_2(1:imax) + ... error correction [2]
                        abs((1:imax) - imax)' / (imax-1) * mltpl + 0.1;
 sec_i = round((length(Tv_exp) - imax) / 4 * 3);
 Tv_exp_2(imax+1:imax+sec_i) = Tv_exp_2(imax+1:imax+sec_i) + ...
                            abs((1:sec_i) - sec_i)' / (sec_i - 1) * 0.1;
 p_e1 = errorbar(shatalov(ind_ini).Tv(:, 1), ... exp plot
   shatalov(ind_ini).Tv(:, 2), Tv_exp_2 .* shatalov(ind_ini).Tv(:, 2), ...
     'ksq', 'MarkerFaceColor', 'k', 'CapSize', 4.5, 'LineWidth', 0.6, ...
                                'DisplayName', "Ibraguimova et al. 2013");
    % plotting the benchmark [3]
 model_vibr_bm  = 2; % benchmark vibrational model
 model_U_bm     = 4; % benchmark U parameter
 p_bm = plot(out(ind_ini, model_vibr_bm, model_U_bm).time_mus, ...
   out(ind_ini, model_vibr_bm, model_U_bm).Tv, pl_vibr(model_vibr_bm), ...
   'Color', pl_U(model_U_bm, :), 'LineWidth', 1.5, 'DisplayName', ...
  "Benchmark: " + str_vibr(model_vibr_bm) + ", U = " + str_U(model_U_bm));
    % plotting simulation results
  for model_U   = 4 %2:4     % choosing U parameter for plotting
 for model_vibr = 2 % 1:2%3  % choosing vibrational models to plot
      % the plot
   p3 = plot(out(ind_ini, model_vibr, model_U).time_mus, ...
       out(ind_ini, model_vibr, model_U).Tv, pl_vibr(model_vibr), ...
       'Color', pl_U(model_U, :), 'LineWidth', 2, ...
         'DisplayName', str_vibr(model_vibr) + ", U = " + str_U(model_U));
      % some plot's parameters
   set(gca, 'FontSize', FontSize);
   ylabel('{\itT}_v, K')%,'FontSize',FontSize);
   xlabel('time, \musec')%,'FontSize',FontSize);
   xlim([min(shatalov(ind_ini).Tv(:, 1)) ... x scale limits
                                max(shatalov(ind_ini).Tv(:, 1)) + 0.2]);
   set(gca, 'TickDir', 'In', ... ticks inside the frame
        'LineWidth', 1, ... width of the frame
        'GridAlpha', 0.2) % grid transparency 
  end
 end
 title("T = " + fix(out(ind_ini, model_vibr, model_U).T(1)) + " K")
 legend('Location', 'best');
end