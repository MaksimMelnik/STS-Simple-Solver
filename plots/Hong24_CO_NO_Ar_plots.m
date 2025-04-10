save_flag=false;

SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output.mat').dat1;

for experiment=[1 2 3]
    mixture_name = sprintf('Mixture%d.csv', experiment);
    path = sprintf('../data/CO_N2_Ar Hong experiment/%s', mixture_name);
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(path);  % Reads CSV into a table
    cleanData = rmmissing(data);        % Remove rows with any NaN

    figure;
    hold on;
    exp_time = cleanData.TimeMs;
    exp_TvCO = cleanData.TvibMeanK;
    plot(exp_time, exp_TvCO, 'linewidth', 1.5, 'DisplayName', 'Экспериментальные данные');

    Tv_CO_SSH = SimulatedData(1,experiment,2).TvCO;
    time_ms_SSH = SimulatedData(1,experiment,2).time/1e3;
    plot(time_ms_SSH, Tv_CO_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');

    Tv_CO_FHO = SimulatedData(2,experiment,2).TvCO;
    time_ms_FHO = SimulatedData(2,experiment,2).time/1e3;
    plot(time_ms_FHO, Tv_CO_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');

    xlabel('Время / мс');
    ylabel('T, K');
    title(['CO-T_{vib}', ' для смеси ', num2str(experiment)]);
    legend('Location', 'southeast');

    lbound = 0.02;
    max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
        time_ms_FHO(time_ms_FHO>lbound), Tv_CO_FHO((time_ms_FHO>lbound)), ...
        'FHO', experiment)
    max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
        time_ms_SSH(time_ms_SSH>lbound), Tv_CO_SSH((time_ms_SSH>lbound)), ...
        'SSH', experiment)
    if save_flag
    end
end

function max_err(time1, values1, time2, values2, model, experiment)
% Max error calculation
% Sample data

% Create common time axis (union of both time vectors)
common_time = union(time1, time2);

% Interpolate both datasets to common time points
interp1_values1 = interp1(time1, values1, common_time, 'linear', 'extrap');
interp1_values2 = interp1(time2, values2, common_time, 'linear', 'extrap');

% Calculate absolute errors
abs_errors = abs(interp1_values1 - interp1_values2);

% Find maximum error
[max_error, max_idx] = max(abs_errors);
fprintf('Experiment %d %s maximum abs error is %.4f at time %.2f\n', experiment, model, ...
    max_error, common_time(max_idx));

rel_errors = abs(interp1_values2 ./ interp1_values1);
% Find maximum error
[max_error, max_idx] = max((1 - rel_errors)*100);
fprintf('Experiment %d %s maximum rel error is %.2f%% at time %.2f\n', experiment, model, ...
    max_error, common_time(max_idx));
end