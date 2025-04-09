save_flag=false;

SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output.mat').dat1;

for experiment=[1 2 3]
    mixture_name = sprintf('Mixture%d.csv', experiment);
    path = sprintf('../data/CO_N2_Ar Hong experiment/%s', mixture_name);
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

    if save_flag
    end
end