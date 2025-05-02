save_flag=false;
drawFlag=false;

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
                model_name = 'Верификационный метод';
        end
        fprintf('\t%s\n', model_name);
        mixture_name = sprintf('Mixture%d.csv', experiment);
        path = sprintf('../data/CO_N2_Ar Hong experiment/%s', mixture_name);
        warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
        data = readtable(path);  % Reads CSV into a table
        cleanData = rmmissing(data);        % Remove rows with any NaN

        
        exp_time = cleanData.TimeMs;
        exp_TvCO = cleanData.TvibMeanK;
        err_index = find(abs(exp_time - 0.1) < 1e-3, 1);
        

        Tv_CO_SSH = SimulatedData(1,experiment, model).TvCO;
        time_ms_SSH = SimulatedData(1,experiment, model).time/1e3;
        

        Tv_CO_FHO = SimulatedData(2,experiment, model).TvCO;
        time_ms_FHO = SimulatedData(2,experiment, model).time/1e3;

        if drawFlag == true
            figure;
            hold on;
            plot(exp_time, exp_TvCO, 'linewidth', 1.5, ...
                'DisplayName', 'Экспериментальные данные');
            errorbar(0.1, exp_TvCO(err_index), err(experiment), 'linewidth', 1.5, ...
                'Color', "#0072BD", ...
                'DisplayName', 'Ошибка в эксперименте');
            plot(time_ms_SSH, Tv_CO_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');
            plot(time_ms_FHO, Tv_CO_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');
            xlabel('Время / мс');
            ylabel('CO-T_{vib}, K');
            %title(['CO-T_{vib}', ' для смеси ', num2str(experiment), '. ', model_name]);
            title(model_name);
            legend('Location', 'southeast');
        end

        if save_flag == true
            filename = sprintf('CO_NO_Ar_Exp%d_Model%d.png', experiment, model);
            exportgraphics(gcf, filename, 'Resolution', 300);
        end

        lbound = 0.06;
        max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
            time_ms_FHO(time_ms_FHO>lbound), Tv_CO_FHO((time_ms_FHO>lbound)), ...
            'FHO', experiment)
        max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
            time_ms_SSH(time_ms_SSH>lbound), Tv_CO_SSH((time_ms_SSH>lbound)), ...
            'SSH', experiment)

        threshold = 0.95;
        print_vibr_relax_time(threshold, exp_time, exp_TvCO, Tv_CO_FHO, time_ms_FHO, 'FHO')
        print_vibr_relax_time(threshold, exp_time, exp_TvCO, Tv_CO_SSH, time_ms_SSH, 'SSH')
    end
end

function print_vibr_relax_time(threshold_const, exp_time, exp, simdata, simtime, model)
    maxT = max(exp);
    lower_bound = threshold_const * maxT;
    valid_T = exp(exp > lower_bound & exp <= maxT);
    average_T = mean(valid_T);
    exp_idx = find(exp >= average_T, 1, 'first');
    fprintf('Exp vibr relax time = %.2f\n', exp_time(exp_idx));
    idx = find(simdata >= average_T, 1, 'first');  % Index of first occurrence
    if isempty(idx)
        fprintf('%s vibr relax time > %.2f\n', model, max(simtime))
    else
        vibr_t = simtime(idx);  % Time corresponding to that index
        fprintf('%s vibr relax time = %.2f\n', model, vibr_t)
        tv_exp = exp_time(exp_idx);
        fprintf('%s vibr relax time error = %.2f%%\n', model, 100*abs((vibr_t-tv_exp)/tv_exp))
    end
end

function max_err(time1, values1, time2, values2, model, experiment)
% Max error calculation
% Sample data

% Create common time axis (union of both time vectors)
common_time = union(time1, time2);

% Interpolate both datasets to common time points
interp1_values1 = interp1(time1, values1, common_time, 'linear', 'extrap');
if experiment == 2 || experiment == 3 % fix for [0.4 0.5] time interval
    interp1_values1(end-4:end) = interp1_values1(end-9:end-5);
end
interp1_values2 = interp1(time2, values2, common_time, 'linear', 'extrap');

% Calculate absolute errors
% abs_errors = abs(interp1_values1 - interp1_values2);

% Find maximum abs error
% [max_error, max_idx] = max(abs_errors);
% fprintf('%s maximum abs error is %.1f at time %.2f\n', model, ...
%     max_error, common_time(max_idx));

rel_errors = abs(interp1_values2 - interp1_values1) ./ interp1_values1;
% Find maximum rel error
[max_error, max_idx] = max(rel_errors*100);
fprintf('%s maximum rel error is %.1f%% at time %.2f\n', model, ...
    max_error, common_time(max_idx));
end