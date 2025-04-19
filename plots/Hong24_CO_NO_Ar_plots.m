save_flag=false;

SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output.mat').dat1;

err = [75 82 64];

for experiment=[1 2 3]
    fprintf('\n\t\tЭксперимент %d\n', experiment);
    for model = [1 2 3]
        switch model
            case 1
                model_name = 'Модель замороженной релаксации';
            case 2
                model_name = 'Модель частичной релаксации';
            case 3
                model_name = 'Верификационная модель';
        end
        fprintf('\t%s\n', model_name);
        mixture_name = sprintf('Mixture%d.csv', experiment);
        path = sprintf('../data/CO_N2_Ar Hong experiment/%s', mixture_name);
        warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
        data = readtable(path);  % Reads CSV into a table
        cleanData = rmmissing(data);        % Remove rows with any NaN

        figure;
        hold on;
        exp_time = cleanData.TimeMs;
        exp_TvCO = cleanData.TvibMeanK;
        err_index = find(abs(exp_time - 0.1) < 1e-3, 1);

        plot(exp_time, exp_TvCO, 'linewidth', 1.5, ...
            'DisplayName', 'Экспериментальные данные');
        errorbar(0.1, exp_TvCO(err_index), err(experiment), 'linewidth', 1.5, ...
            'Color', "#0072BD", ...
            'DisplayName', 'Ошибка в эксперименте');

        Tv_CO_SSH = SimulatedData(1,experiment, model).TvCO;
        time_ms_SSH = SimulatedData(1,experiment, model).time/1e3;
        plot(time_ms_SSH, Tv_CO_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');

        Tv_CO_FHO = SimulatedData(2,experiment, model).TvCO;
        time_ms_FHO = SimulatedData(2,experiment, model).time/1e3;
        plot(time_ms_FHO, Tv_CO_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');

        xlabel('Время / мс');
        ylabel('T, K');
        title(['CO-T_{vib}', ' для смеси ', num2str(experiment), '. ', model_name]);
        legend('Location', 'southeast');
        if save_flag == true
            filename = sprintf('CO_NO_Ar_Exp%d_Model%d.png', experiment, model);
            exportgraphics(gcf, filename, 'Resolution', 300);
        end

        lbound = 0.01;
        max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
            time_ms_FHO(time_ms_FHO>lbound), Tv_CO_FHO((time_ms_FHO>lbound)), ...
            'FHO')
        max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
            time_ms_SSH(time_ms_SSH>lbound), Tv_CO_SSH((time_ms_SSH>lbound)), ...
            'SSH')
    end
end

function max_err(time1, values1, time2, values2, model)
% Max error calculation
% Sample data

% Create common time axis (union of both time vectors)
common_time = union(time1, time2);

% Interpolate both datasets to common time points
interp1_values1 = interp1(time1, values1, common_time, 'linear', 'extrap');
interp1_values2 = interp1(time2, values2, common_time, 'linear', 'extrap');

% Calculate absolute errors
% abs_errors = abs(interp1_values1 - interp1_values2);

% Find maximum abs error
% [max_error, max_idx] = max(abs_errors);
% fprintf('%s maximum abs error is %.1f at time %.2f\n', model, ...
%     max_error, common_time(max_idx));

rel_errors = abs(interp1_values2 ./ interp1_values1);
% Find maximum rel error
[max_error, max_idx] = max((1 - rel_errors)*100);
fprintf('%s maximum rel error is %.1f%% at time %.2f\n', model, ...
    max_error, common_time(max_idx));
end