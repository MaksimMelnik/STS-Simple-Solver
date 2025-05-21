save_flag=true;
drawFlag=false;

% errors dictionary
err = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
err(8) = 75;
err(10) = 82;
err(4) = 64;

vt_only = true;
dbg = false;

if vt_only
    SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output_VT.mat').result;
else
    SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output_VT_VV.mat').result;
end

for mixture=1:3
    fprintf('\n\t\t\tСмесь %d\n', mixture);

    max_exp = 10;
    if mixture == 3
        max_exp = 14;
    end

    for experiment=1:max_exp
        fprintf('\t\tЭксперимент %d\n', experiment);
        for model = 1:3
            switch model
                case 1
                    model_name = 'Модель замороженной релаксации';
                case 2
                    model_name = 'Модель частичной релаксации';
                case 3
                    model_name = 'Верификационный метод';
            end
            fprintf('\t%s\n', model_name);
            path = sprintf('../data/CO_N2_Ar Hong experiment/Mixture%d/Experiment%d', mixture, experiment);

            opts = detectImportOptions(path);
            opts.VariableNamesLine = 11;
            opts.DataLines = [12, Inf];
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
            data = readtable(path, opts);  % Reads CSV into a table
            cleanData = rmmissing(data, 'MinNumMissing', 6);        % Remove rows with any NaN


            exp_time = cleanData.TimeMs;
            exp_TvCO = cleanData.TvibMeanK;
            err_index = find(abs(exp_time - 0.1) < 1e-3, 1);
            Tv_CO_SSH = SimulatedData(mixture,experiment,model,1).TvCO;
            time_ms_SSH = SimulatedData(mixture,experiment,model,1).time/1e3;
            Tv_CO_FHO = SimulatedData(mixture,experiment,model,2).TvCO;
            time_ms_FHO = SimulatedData(mixture,experiment,model,2).time/1e3;

            if drawFlag
                figure;
            else
                figTvCO = figure('Visible', 'off');
            end
            hold on;
            plot(exp_time, exp_TvCO, 'linewidth', 1.5, ...
                'DisplayName', 'Экспериментальные данные');
            if calc_errors(mixture, experiment)
                errorbar(0.1, exp_TvCO(err_index), err(experiment), 'linewidth', 1.5, ...
                    'Color', "#0072BD", ...
                    'DisplayName', 'Ошибка в эксперименте');
            end
            plot(time_ms_SSH, Tv_CO_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');
            plot(time_ms_FHO, Tv_CO_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');
            xlabel('Время / мс');
            ylabel('CO-T_{vibr}, K');
            title(sprintf('Смесь %d эксперимент %d: %s', mixture, experiment, model_name));
            legend('Location', 'southeast');


            if save_flag == true
                filename = sprintf('../data/CO_N2_Ar Hong experiment/Mixture%d/Tv_CO_Exp%d_Model%d.png', mixture, experiment, model);
                exportgraphics(figTvCO, filename, 'Resolution', 300);
            end

            exp_time_2 = cleanData.TimeMs_2(cleanData.TimeMs_2 > 0);
            exp_p = cleanData.PressureBar(cleanData.TimeMs_2 > 0);
            p_SSH = SimulatedData(mixture,experiment,model,1).p * 133.322368 / 1e5; % Pa to Bar
            p_FHO = SimulatedData(mixture,experiment,model,2).p * 133.322368 / 1e5;
            if drawFlag
                figure;
            else
                figPressureP = figure('Visible', 'off');
            end
            hold on;
            plot(exp_time_2, exp_p, 'linewidth', 1.5, ...
                'DisplayName', 'Экспериментальные данные');
            plot(time_ms_SSH, p_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');
            plot(time_ms_FHO, p_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');
            xlabel('Время / мс');
            ylabel('p, бар');
            title(sprintf('Смесь %d эксперимент %d: %s', mixture, experiment, model_name));
            legend('Location', 'southeast');

            if save_flag == true
                filename = sprintf('../data/CO_N2_Ar Hong experiment/Mixture%d/P_Exp%d_Model%d.png', mixture, experiment, model);
                exportgraphics(figPressureP, filename, 'Resolution', 300);
            end

            exp_time_1 = cleanData.TimeMs(cleanData.TimeMs > 0);
            exp_T = cleanData.TrotMeanK(cleanData.TimeMs > 0);
            T_SSH = SimulatedData(mixture,experiment,model,1).T;
            T_FHO = SimulatedData(mixture,experiment,model,2).T;
            if drawFlag
                figure;
            else
                figPressureT = figure('Visible', 'off');
            end
            hold on;
            plot(exp_time_1, exp_T, 'linewidth', 1.5, ...
                'DisplayName', 'Экспериментальные данные');
            plot(time_ms_SSH, T_SSH, 'g:', 'linewidth', 2, 'DisplayName', 'SSH');
            plot(time_ms_FHO, T_FHO, 'r--', 'linewidth', 2, 'DisplayName', 'FHO');
            xlabel('Время / мс');
            ylabel('T, К');
            title(sprintf('Смесь %d эксперимент %d: %s', mixture, experiment, model_name));
            legend('Location', 'southeast');

            if save_flag == true
                filename = sprintf('../data/CO_N2_Ar Hong experiment/Mixture%d/T_Exp%d_Model%d.png', mixture, experiment, model);
                exportgraphics(figPressureT, filename, 'Resolution', 300);
            end

            if calc_errors(mixture, experiment)
                lbound = 0.01;
                max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
                    time_ms_FHO(time_ms_FHO>lbound), Tv_CO_FHO((time_ms_FHO>lbound)), ...
                    'FHO')
                max_err(exp_time(exp_time>lbound), exp_TvCO(exp_time>lbound), ...
                    time_ms_SSH(time_ms_SSH>lbound), Tv_CO_SSH((time_ms_SSH>lbound)), ...
                    'SSH')

                threshold = 0.95;
                print_vibr_relax_time(threshold, exp_time, exp_TvCO, Tv_CO_FHO, time_ms_FHO, 'FHO', experiment)
                print_vibr_relax_time(threshold, exp_time, exp_TvCO, Tv_CO_SSH, time_ms_SSH, 'SSH', experiment)
            end
        end
    end
end

function [flag] = calc_errors(mixture, experiment)
flag = false;
if mixture == 1 && experiment == 8
    flag = true;
end
if mixture == 2 && experiment == 10
    flag = true;
end
if mixture == 3 && experiment == 4
    flag = true;
end
end

function print_vibr_relax_time(threshold_const, exp_time, exp_T, sim_T, sim_time, model, exp1_idx)
global dbg;
maxT = max(exp_T);
lower_bound = threshold_const * maxT;
T_valid = exp_T(exp_T > lower_bound & exp_T <= maxT);
T_v_eq = mean(T_valid);
exp_idx = find(exp_T >= T_v_eq, 1, 'first');
idx = find(sim_T >= T_v_eq, 1, 'first');  % Index of first occurrence

if dbg == true
    fprintf('Exp full equilibrium relaxation time = %.2f\n', exp_time(exp_idx));
    if isempty(idx)
        fprintf('%s Full equilibrium time > %.2f\n', model, max(sim_time))
    else
        vibr_t = sim_time(idx);  % Time corresponding to that index
        fprintf('%s Full equilibrium time = %.2f\n', model, vibr_t)
        tv_exp = exp_time(exp_idx);
        fprintf('%s Full equilibrium time error = %.2f%%\n', model, 100*abs((vibr_t-tv_exp)/tv_exp))
    end
end

switch exp1_idx
    case 8
        correct_Tv_0_idx = 83;
    case 10
        correct_Tv_0_idx = 79;
    case 4
        correct_Tv_0_idx = 88;
    otherwise
        error('Bad index');
end
T_0_exp = exp_T(correct_Tv_0_idx);
T_relax_exp = (1/exp(1)) * T_0_exp + (1 - 1/exp(1)) * T_v_eq;
[~, idx_tau_exp] = min(abs(exp_T - T_relax_exp));
tau_exp = exp_time(idx_tau_exp);
fprintf('Exp tau = %.3f\n', tau_exp);

maxT = max(sim_T);
lower_bound = threshold_const * maxT;
T_valid = sim_T(sim_T > lower_bound & sim_T <= maxT);
T_v_eq = mean(T_valid);

T_0_model = sim_T(1);
T_relax_model = (1/exp(1)) * T_0_model + (1 - 1/exp(1)) * T_v_eq;
[~, idx_tau_model] = min(abs(sim_T - T_relax_model));
tau_model = sim_time(idx_tau_model);
fprintf('%s tau = %.3f\n', model, tau_model);
fprintf('%s tau error = %.1f%%\n', model, 100 * abs(tau_exp - tau_model)/tau_exp);
end

function max_err(time1, values1, time2, values2, model)
% Max error calculation
% Sample data

% Create common time axis (union of both time vectors)
values2 = values2(time2 < max(time1));
time2 = time2(time2 < max(time1));

common_time = union(time1, time2);

% Interpolate both datasets to common time points
interp1_values1 = interp1(time1, values1, common_time, 'linear', 'extrap');
[max_val, max_idx] = max(interp1_values1);
interp1_values1(max_idx+1:end) = max_val;
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