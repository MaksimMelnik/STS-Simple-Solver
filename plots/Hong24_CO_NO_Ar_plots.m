save_flag=false;

SimulatedData = load('../data/CO_N2_Ar Hong experiment/CO_N2_Ar_behindRSW_output.mat').dat1;

for experiment=[1 2 3]
    mixture_name = sprintf('Mixture%d.csv', experiment);
    path = sprintf('../data/CO_N2_Ar Hong experiment/%s', mixture_name);
    data = readtable(path);  % Reads CSV into a table
    cleanData = rmmissing(data);        % Remove rows with any NaN
    exp_time = cleanData.TimeMs;
    exp_TvCO = cleanData.TvibMeanK;

    for kinetic_model = [2] % [SSH, FHO]
        Tv_CO = SimulatedData(kinetic_model,experiment,2).TvCO;
        time_ms = SimulatedData(kinetic_model,experiment,2).time/1e3;
        figure

        plot(time_ms, Tv_CO, exp_time, exp_TvCO, 'linewidth', 1.5);
        xlabel('Время / мс');
        ylabel('T, K');
        title(['CO-T_{vib}', ' for Mixture ', num2str(experiment)]);
        if save_flag
        end
    end
end