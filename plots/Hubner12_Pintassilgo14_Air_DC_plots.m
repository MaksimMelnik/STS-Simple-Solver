function Hubner12_Pintassilgo14_Air_DC_plots(data)
% Comparison plots for simulations of Hubner experiment [1] with previous
% simulations by Pintassilgo [2].
% [1] M Hubner et al Meas. Sci. Technol. 23 (2012) 115602.
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.
% 26.07.2023 Maksim Melnik

is_T_Tv_plot = true;
is_Q_plot = true;
is_n_NO_N_O = false;
is_VDF = true;
is_Q_plot2 = true;
is_Omega_plot = false;
is_Pintassilgo2014_T = true;

    % constants
N_a = 6.02214076e23;              % Avogadro constant
k = 1.380649e-23;                 % Boltzmann constant, J/K
    % loading data from papers
load('../data/for comparison/Hubner2012_and_Pintassilgo2014.mat' ...
                                                            ) %#ok<LOAD>
addpath('../src/')

    % taking required data
kinetics = data.kinetics;
Y  = data.Y;
t  = data.t;
T  = data.T;
if isKey(kinetics.reactions, 'free_e')
    Te = data.Te;
    Free_e_reactions = kinetics.reactions("free_e");
    IndexOfMolecules = kinetics.IndexOfMolecules;
    e.name = "e";
    e.form_e = 0;
    e.num_vibr_levels = 0;
    e.e_E             = 0;
    e.fr_deg_c        = 3;
    e.s_e             = 1;
    IndexOfMolecules("e") = kinetics.num_Ps + 1;
    index                 = [kinetics.index, kinetics.num_eq + 1];
    Ps                    = [kinetics.Ps, e];
end
Tv = data.Tv;
N2 = kinetics.Ps{1};
O2 = kinetics.Ps{2};
NO = kinetics.Ps{3};
N  = kinetics.Ps{4};
O  = kinetics.Ps{5};
if isKey(kinetics.reactions, 'Exch')
    Exch_reactions = kinetics.reactions("Exch");
    IndexOfMolecules = kinetics.IndexOfMolecules;
end

n_g=sum(Y(:, 1:end-1), 2);
n = {};
for i = 1:length(kinetics.Ps)
    n{i} = sum(Y(:, kinetics.index{i}), 2);
end
n_N2=n{1};
t_ag=t-0.005;
fsize = [200 50 900 550];
Exch_colors = [
                41, 128, 185
                192, 57, 43
                39, 174, 96
                243, 156, 18
                241, 196, 15
                142, 68, 173
                26, 188, 156
                234, 76, 137
                65, 0, 147
                160, 82, 45
                50, 205, 50
                ]/255;

%%  Omega
if is_Omega_plot

%VT
for i_M1 = 1:3
    figure
    hold on
    M1 = kinetics.Ps{i_M1};
    iM1 = kinetics.index{i_M1};
    for i_M2 = 1:5
        R_VT_M1 = zeros(length(t), 1);
        M2 = kinetics.Ps{i_M2};
        iM2 = kinetics.index{i_M2};
        for i_out = 1:length(t)
            for ind_e = 1:M1.num_elex_levels
                if M1.num_vibr_levels(ind_e) > 1
                    if ind_e > 1
                        i_e_prev = M1.num_vibr_levels(ind_e - 1);
                    else
                        i_e_prev = 0;
                    end
                    ie_M1 = iM1(i_e_prev + 1:i_e_prev + M1.num_vibr_levels(ind_e));
                    [R_VT_data, ~] = R_VT(M1, Y(i_out, ie_M1)', M2, ...
                        sum(Y(i_out, iM2)), T(i_out), ind_e, kinetics.reactions('VT'));
                    R_VT_M1(i_out) = R_VT_M1(i_out) + sum(R_VT_data .* (M1.ev_0(ind_e) + M1.ev_i{ind_e} + M1.form_e)');
                end
            end
        end
        plot(t*1e3, R_VT_M1, 'linewidth', 1.5);
    end
    s = strcat('\Omega VT, ', M1.name());
    title(s);
    legend('+ N2', '+ O2', '+ NO', '+ N', '+ O');
    xlim([0 20]);
    hold off
end

%VV
for i_M1 = 1:3
    figure
    hold on
    M1 = kinetics.Ps{i_M1};
    iM1 = kinetics.index{i_M1};
    for i_M2 = 1:3
        R_VV_M1 = zeros(length(t), 1);
        M2 = kinetics.Ps{i_M2};
        iM2 = kinetics.index{i_M2};
        for i_out = 1:length(t)
            for ind_e = 1:M1.num_elex_levels
                if M1.num_vibr_levels(ind_e) > 1
                    if ind_e > 1
                        i_e_prev = M1.num_vibr_levels(ind_e - 1);
                    else
                        i_e_prev = 0;
                    end
                    ie_M1 = iM1(i_e_prev + 1:i_e_prev + M1.num_vibr_levels(ind_e));
                    for ind_eM2 = 1:M2.num_elex_levels
                        if M2.num_vibr_levels(ind_eM2) > 1
                            if ind_eM2 > 1
                                i_e_prev = M2.num_vibr_levels(ind_eM2 - 1);
                            else
                                i_e_prev = 0;
                            end
                            ie_M2 = iM2(i_e_prev + 1:i_e_prev + M2.num_vibr_levels(ind_eM2));
                            [R_VV_data, ~] = R_VV(M1, Y(i_out, ie_M1)', M2, Y(i_out, ie_M2)', ...
                                    T(i_out), ind_e, ind_eM2, kinetics.reactions('VV'));
                            R_VV_M1(i_out) = R_VV_M1(i_out) + sum((M1.ev_0(ind_e) + M1.ev_i{ind_e} + M1.form_e) * R_VV_data);
                        end
                    end
                end
            end
        end
        plot(t*1e3, R_VV_M1, 'linewidth', 1.5); %
    end
    s = strcat('\Omega VV, ', M1.name());
    title(s);
    legend('+ N2', '+ O2', '+ NO');
    xlim([0 20]);
    hold off
end

%Exch
R_exch_data_full = zeros(length(kinetics.Ps), length(Exch_reactions), length(t));
for ind_exch = 1:length(Exch_reactions)
    reaction = Exch_reactions(ind_exch);
    
    IOM_M1 = IndexOfMolecules(reaction.particles(1));
    IOM_M2 = IndexOfMolecules(reaction.particles(2));
    IOM_M3 = IndexOfMolecules(reaction.particles(3));
    IOM_M4 = IndexOfMolecules(reaction.particles(4));
    M1 = kinetics.Ps{IOM_M1};
    M2 = kinetics.Ps{IOM_M2};
    M3 = kinetics.Ps{IOM_M3};
    M4 = kinetics.Ps{IOM_M4};
    indM1 = kinetics.index{IOM_M1};
    indM2 = kinetics.index{IOM_M2};
    indM3 = kinetics.index{IOM_M3};
    indM4 = kinetics.index{IOM_M4};
    ie_1 = reaction.index{1}{1};
    ie_2 = reaction.index{2}{1};
    ie_3 = reaction.index{3}{1};
    ie_4 = reaction.index{4}{1};
    for i_out = 1:length(t)
     switch length(reaction.particles)
      case 4
       [R_exch_temp, ~] = R_exch({M1, M2, M3, M4}, Y(i_out, indM1)', ...
           Y(i_out, indM2)', Y(i_out, indM3)', Y(i_out, indM4)', ...
           T(i_out), reaction);
       if M1.num_vibr_levels(ie_1) > 1
            i1_e = (1+sum(M1.num_vibr_levels(1:ie_1-1)) : ...
                                        sum(M1.num_vibr_levels(1:ie_1)));
            R_exch_data_full(IOM_M1, ind_exch, i_out) = ...
                sum(R_exch_temp(i1_e, :, :, :) .* (M1.ev_0(ie_1) ...
                + M1.ev_i{ie_1}' + M1.e_E(ie_1) + M1.form_e), 'all');
       end
       if M2.num_vibr_levels(ie_2) > 1
            i2_e = (1+sum(M2.num_vibr_levels(1:ie_2 - 1)) : ...
                                        sum(M2.num_vibr_levels(1:ie_2)));
            R_exch_data_full(IOM_M2, ind_exch, i_out) = ...
                sum(R_exch_temp(:, i2_e, :, :) .* (M2.ev_0(ie_2) ...
                + M2.ev_i{ie_2} + M2.e_E(ie_2) + M2.form_e), 'all');
       else
            R_exch_data_full(IOM_M2, ind_exch, i_out) = ...
                sum(R_exch_temp(:, 1, :, :) * ...
                (M2.e_E(ie_2) + M2.form_e), 'all');
       end
       if M3.num_vibr_levels(ie_3) > 1
            i3_e = (1+sum(M3.num_vibr_levels(1:ie_3 - 1)) : ...
                                        sum(M3.num_vibr_levels(1:ie_3)));
            R_exch_data_full(IOM_M3, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, i3_e , :) .* ...
                reshape(M3.ev_0(ie_3) + M3.ev_i{ie_3}' + M3.e_E(ie_3) ...
                + M3.form_e, 1, 1, []), 'all');
       end
       if M4.num_vibr_levels(ie_4) > 1
            i4_e = (1+sum(M4.num_vibr_levels(1:ie_4 - 1)) : ...
                                        sum(M4.num_vibr_levels(1:ie_4)));
            R_exch_data_full(IOM_M4, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, i4_e) .* ...
                reshape(M4.ev_0(ie_4) + M4.ev_i{ie_4}' + M4.e_E(ie_4) ...
                + M4.form_e, 1, 1, 1, []), 'all');
       else
            R_exch_data_full(IOM_M4, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, 1) * ...
                (M4.e_E(ie_4) + M4.form_e), 'all');
       end
      case 5
       IOM_M5 = IndexOfMolecules(reaction.particles(5));
       indM5  = kinetics.index{IOM_M5};
       M5 = kinetics.Ps{IOM_M5};
       ie_5 = reaction.index{5}{1};
       [R_exch_temp, ~] = R_exch_23({M1, M2, M3, M4, M5}, ...
            Y(i_out, indM1)', Y(i_out, indM2)', Y(i_out, indM3)', ...
            Y(i_out, indM4)', Y(i_out, indM5)', T(i_out), reaction, 1);
        if M1.num_vibr_levels(ie_1) > 1
            i1_e = (1+sum(M1.num_vibr_levels(1:ie_1 - 1)) : ...
                                        sum(M1.num_vibr_levels(1:ie_1)));
            R_exch_data_full(IOM_M1, ind_exch, i_out) = ...
                sum(R_exch_temp(i1_e, :, :, :, :) .* (M1.ev_0(ie_1) ...
                + M1.ev_i{ie_1}' + M1.e_E(ie_1) + M1.form_e), 'all');
        end
        if M2.num_vibr_levels(ie_2) > 1
            i2_e = (1+sum(M2.num_vibr_levels(1:ie_2 - 1)) : ...
                                        sum(M2.num_vibr_levels(1:ie_2)));
            R_exch_data_full(IOM_M2, ind_exch, i_out) = ...
                sum(R_exch_temp(:, i2_e , :, :, :) .* (M2.ev_0(ie_2) ...
                + M2.ev_i{ie_2} + M2.e_E(ie_2) + M2.form_e), 'all');
        else
            R_exch_data_full(IOM_M2, ind_exch, i_out) = ...
                sum(R_exch_temp(:, 1, :, :, :) * ...
                (M2.e_E(ie_2) + M2.form_e), 'all');
        end
        if M3.num_vibr_levels(ie_3) > 1
            i3_e = (1+sum(M3.num_vibr_levels(1:ie_3 - 1)) : ...
                                        sum(M3.num_vibr_levels(1:ie_3)));
            R_exch_data_full(IOM_M3, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, i3_e, :, :) .* ...
                reshape(M3.ev_0(ie_3) + M3.ev_i{ie_3}' + M3.e_E(ie_3) ...
                + M3.form_e, 1, 1, []), 'all');
        end
        if M4.num_vibr_levels(ie_4) > 1
            i4_e = (1+sum(M4.num_vibr_levels(1:ie_4 - 1)) : ...
                                        sum(M4.num_vibr_levels(1:ie_4)));
            R_exch_data_full(IOM_M4, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, i4_e, :) .*...
                reshape(M4.ev_0(ie_4) + M4.ev_i{ie_4}' + M4.e_E(ie_4) ...
                + M4.form_e, 1, 1, 1, []), 'all');
        else
            R_exch_data_full(IOM_M4, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, 1, :) * ...
                (M4.e_E(ie_4) + M4.form_e), 'all');
        end
        if M5.num_vibr_levels(ie_5) > 1
            i5_e = (1+sum(M5.num_vibr_levels(1:ie_5 - 1)) : ...
                                        sum(M5.num_vibr_levels(1:ie_5)));
            R_exch_data_full(IOM_M5, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, :, i5_e) .* ...
                reshape(M5.ev_0(ie_5) + M5.ev_i{ie_5}' + M5.e_E(ie_5) ...
                + M5.form_e, 1, 1, 1, 1, []), 'all');
        else
            R_exch_data_full(IOM_M5, ind_exch, i_out) = ...
                sum(R_exch_temp(:, :, :, :, 1) * ...
                (M5.e_E(ie_5) + M5.form_e), 'all');
        end
      otherwise
            error("The number of particles is not 4 or 5.")
     end
    end
end

for ind = 1:5
figure
hold on
l = {};
for ind_exch = 1:length(Exch_reactions)
    plot(t*1e3, reshape(R_exch_data_full(ind, ind_exch, :), [], 1), ...
                    'linewidth', 1.5, 'Color', Exch_colors(ind_exch, :));
    l{ind_exch} = Exch_reactions(ind_exch).name; 
end
title(strcat("\Omega exchange reactions, ", kinetics.Ps{ind}.name));
xlim([0 20]);
legend(l);
hold off
end

% Free e
R_free_e_data_full = zeros(length(kinetics.Ps), length(Free_e_reactions), length(t));
for ind_free_e = 1:length(Free_e_reactions)
    reaction = Free_e_reactions(ind_free_e);

    IOM_M1 = IndexOfMolecules(reaction.particles(1));
    IOM_M2 = IndexOfMolecules(reaction.particles(2));
    IOM_M3 = IndexOfMolecules(reaction.particles(3));
    IOM_M4 = IndexOfMolecules(reaction.particles(4));
    M1 = Ps{IOM_M1};
    M2 = Ps{IOM_M2};
    M3 = Ps{IOM_M3};
    M4 = Ps{IOM_M4};
    indM1 = index{IOM_M1};
    indM2 = index{IOM_M2};
    indM3 = index{IOM_M3};
    indM4 = index{IOM_M4};
    switch length(reaction.particles)
    case 4
        for i_out = 1:length(t)
            [R_free_e_temp, ~] = R_exch({M1, M2, M3, M4}, ...
                        Y(i_out, indM1)', Y(i_out, indM2)', Y(i_out, indM3)', Y(i_out, indM4)', T(i_out), reaction);
            R_free_e_data_full(IOM_M1, ind_free_e, i_out) = R_free_e_temp * (3/2*k*Te(i_out));
            R_free_e_data_full(IOM_M2, ind_free_e, i_out) = R_free_e_temp * (M2.form_e);
            R_free_e_data_full(IOM_M3, ind_free_e, i_out) = -2*R_free_e_temp * (M3.form_e);
        end
    case 5
        IOM_M5 = IndexOfMolecules(reaction.particles(5));
        M5 = Ps{IOM_M5};
        indM5  = kinetics.index{IOM_M5};
        for i_out = 1:length(t)
            [R_free_e_temp, ~] = R_exch_23_e({M1, M2, M3, M4, M5}, ...
                    Y(i_out, indM1)', Y(i_out, indM2)', Y(i_out, indM3)', Y(i_out, indM4)', Y(i_out, indM5)', Te(i_out), reaction);
            R_free_e_data_full(IOM_M2, ind_free_e, i_out) = sum(R_free_e_temp(1, 1:M2.num_vibr_levels(1), :, :) .* (M2.ev_0(1) + M2.ev_i{1} + M2.form_e), 'all');
            R_free_e_data_full(IOM_M4, ind_free_e, i_out) = -2*sum(R_free_e_temp(1, 1:M2.num_vibr_levels(1), :, :)) * (M4.form_e);
        end
    end
    
end

for ind = 1:size(Ps, 2)
    figure
    hold on
    l = {};
    for ind_exch = 1:length(Free_e_reactions)
        plot(t*1e3, reshape(R_free_e_data_full(ind, ind_exch, :), [], 1), 'linewidth', 1.5);
        l{ind_exch} = Free_e_reactions(ind_exch).name;
    end
    title(strcat("\Omega free e reactions, ", Ps{ind}.name));
    legend(l);
    xlim([0 20]);
    hold off
end

% elastic collisions e with particles
R_e_data_full = zeros(length(kinetics.Ps), length(t));
figure
hold on
for iM = 1:length(kinetics.Ps)
    M = kinetics.Ps{iM};
    indM = kinetics.index{iM};
    me = 9.1094e-31;
    mi = kinetics.Ps{iM}.mass;
    mred = mi*me/(mi+me);
    r = M.diameter / 2;
    for i_out = 1:length(t)
     ni = sum(Y(i_out, indM));
     ind_e = max(kinetics.index{length(kinetics.Ps)})+1;
     ne = Y(i_out, ind_e);
     z = sqrt(8*pi*k*Te(i_out)/mred)*r^2;                          % m3/s
     R_e_data_full(iM, i_out) = z*ni*ne*me/mi...
                                        * (3/2*k*Te(i_out));
    end
    plot(t*1e3, R_e_data_full(iM, :), 'linewidth', 1.5);
end
title("\Omega elastic collisions e with particles");
legend('+N2', '+O2', '+NO', '+N', '+O', '+N2+', '+O2+');
xlim([0 20]);
hold off
end
%% T and Tv plot
if is_T_Tv_plot 
    figure  
plot(t*1e3, T, t*1e3, Tv, '-.', 'linewidth', 1.5)
legend('T', 'Tv', 'location', 'best')
xlabel('t, ms')
% xlim([6e-2 1.5e-1])
end
%% T vs exp plot
    figure('Position', fsize)
tH = Hubner_2012_T(:, 1);
TH = Hubner_2012_T(:, 2);
hold on
plot(tH, TH, 'sq', 'MarkerEdgeColor', '#5A5A5A', 'linewidth', 1.5);
plot(t*1e3, T, 'linewidth', 1.5);
legend_str = ["T_{exp}, Hubner 2012", "T"];
if is_Pintassilgo2014_T
    plot(Pintassilgo2014_T(:, 1), Pintassilgo2014_T(:, 2), ':', ...
                                                        'LineWidth', 1.5)
    legend_str = [legend_str, "T, Pintassilgo 2014"];
end
err = 30*ones(size(TH(5:5:end)))';
errorbar(tH(5:5:end), TH(5:5:end),...
                     err, 's', 'color', '#5A5A5A', 'linewidth', 1);
hold off
legend(legend_str, 'location', 'best')
xlabel('t, ms')
xlim([-2 14])
ylim([250 620])
%% Te plot
if isKey(kinetics.reactions, 'free_e')
        figure
    plot(t*1e3, Te, 'linewidth', 1.5);
    xlabel('t, ms')
    legend('Te')
end
%% N, O and NO ag plot
    figure('Position', fsize)
loglog(Pintassilgo2014_ag_N(:, 1), Pintassilgo2014_ag_N(:, 2), ...
                        'color', [0.9 0 0], 'linewidth', 1.5) %#ok<USENS>
hold on
loglog(t_ag*1e3, Y(:, kinetics.index{4})./n_g, ...
    ':', 'color', [0.9 0 0], 'linewidth', 1.5)
loglog(Pintassilgo2014_ag_O(:, 1), Pintassilgo2014_ag_O(:, 2), ...
                        'color', [0 0.6 0], 'linewidth', 1.5) %#ok<USENS>
loglog(t_ag*1e3, Y(:, kinetics.index{5})./n_g, ...
    ':', 'color', [0 0.6 0], 'linewidth', 1.5)
loglog(Pintassilgo2014_ag_NO(:, 1), Pintassilgo2014_ag_NO(:, 2), ...
                        'color', [0 0 0.9], 'linewidth', 1.5) %#ok<USENS>
loglog(t_ag*1e3, sum(Y(:, kinetics.index{3}), 2)./n_g, ...
    ':', 'color', [0 0 0.9], 'linewidth', 1.5)
loglog(t_ag*1e3, n_N2./n_g, ...
    ':', 'color', [0 0 0], 'linewidth', 2.5)
legend('n_N, Pintassilgo2014', 'n_N, Maksim', ...
    'n_O, Pintassilgo2014', 'n_O, Maksim', ...
    'n_{NO}, Pintassilgo2014', 'n_{NO}, Maksim', 'N_2, Maksim', ...
    'location', 'best')
xlim([1e-2 2.5e2])
ylim([1e-4 1])
grid on

%% VDF ag plot
if is_VDF
time_ind0=1;
[~, time_ind1]   = min( abs(t_ag*1e3 - 1) );
[~, time_ind10]  = min( abs(t_ag*1e3 - 10) );
% time_ind100=length(Y(:, 1));
[~, time_ind100] = min( abs(t_ag*1e3 - 100) );
for i = 1:length(kinetics.Ps)
if length(Y(1, kinetics.index{i}))-1 ~= 0
    lvls4plot=0:length(Y(1, kinetics.index{i}))-1;
    figure('Position', fsize)
    semilogy(lvls4plot, Y(time_ind0, kinetics.index{i})/n{i}(time_ind0), ...
                            '--', 'color', [1 0 0] ,'linewidth', 2.5)
    hold on
    semilogy(lvls4plot, Y(time_ind1, kinetics.index{i})/n{i}(time_ind1), ...
                        '--', 'color', [0.75 0.2 0.25] ,'linewidth', 2.5)
    semilogy(lvls4plot, Y(time_ind10, kinetics.index{i})/n{i}(time_ind10), ...
                        '--', 'color', [0.25 0.4 0.75] ,'linewidth', 2.5)
    semilogy(lvls4plot, ...
            Y(time_ind100, kinetics.index{i})/n{i}(time_ind100), ...
                            '--', 'color', [0 0.6 1] ,'linewidth', 2.5)
    if i == 1
        semilogy(Pintassilgo2014_N2_VDF_post_DC(:,1), ...
                        Pintassilgo2014_N2_VDF_post_DC(:,2), ...
                        'color', [1 0 0] ,'linewidth', 1.5) %#ok<USENS>
        semilogy(Pintassilgo2014_N2_VDF_ag_1ms(:,1), ...
                        Pintassilgo2014_N2_VDF_ag_1ms(:,2), ...
                        'color', [0.75 0.2 0.25] ,'linewidth', 1.5) %#ok<USENS>
        semilogy(Pintassilgo2014_N2_VDF_ag_10ms(:,1), ...
                        Pintassilgo2014_N2_VDF_ag_10ms(:,2), ...
                        'color', [0.25 0.4 0.75] ,'linewidth', 1.5) %#ok<USENS>
        semilogy(Pintassilgo2014_N2_VDF_ag_100ms(:,1), ...
                        Pintassilgo2014_N2_VDF_ag_100ms(:,2), ...
                        'color', [0 0.6 1] ,'linewidth', 1.5) %#ok<USENS>
    end
    hold off
    title(kinetics.Ps{i}.name)
    if kinetics.Ps{i}.name == "N2"
        legend('0 ms, Maksim', '1 ms, Maksim', '10 ms, Maksim', ...
            '100 ms, Maksim', '0.1 ms, Pintassilgo2014', ...
            '1 ms, Pintassilgo2014', '10 ms, Pintassilgo2014', ...
            '100 ms, Pintassilgo2014', 'location', 'best')
    else
        legend('0 ms, Maksim', '1 ms, Maksim', '10 ms, Maksim', ...
            '100 ms, Maksim', 'location', 'best')
    end
    xlim([0 30])
    ylim([1e-6 1])
end
end
end
%% heating rates, K/s
if is_Q_plot 
    figure('Position', fsize)
iN2 = kinetics.index{1};
i1_N2 = iN2(1:N2.num_vibr_levels(1));
iO = kinetics.index{5};
Q_VT = zeros(length(t), 1);
% Q_VT_wall = zeros(length(t), 1);
Q_rec_wall = zeros(length(t), 1);
Q_VV = zeros(length(t), 1);
Q_exch_N2_O = zeros(length(t), 1);
Q_exch_N_NO = zeros(length(t), 1);
for i_out = 1:length(t)
if isKey(kinetics.reactions, 'VT')
    [~, Q_VT_data] = R_VT(N2, Y(i_out, i1_N2)', O, ...
                Y(i_out, iO(1)), T(i_out), 1, kinetics.reactions('VT'));
    Q_VT(i_out) = Q_VT_data;
end
if isKey(kinetics.reactions, 'Wall') 
%  [~, Q_VT_wall_data] = ...
%                     R_VT_wall(N2, Y(i_out, i1_N2)', T(i_out), kinetics);
%  Q_VT_wall(i_out) = Q_VT_wall_data;
 if isKey(kinetics.reactions, 'Rec_wall')
  [~, Q_rec_wall_data] = R_rec_wall(O, Y(i_out, kinetics.index{5}), ...
                                                    T(i_out), kinetics);
  Q_rec_wall(i_out) = Q_rec_wall_data;
 end
end
if isKey(kinetics.reactions, 'VV')
 [~, Q_VV_data] = R_VV(N2, Y(i_out, i1_N2)', N2, Y(i_out, i1_N2)', ...
                             T(i_out), 1, 1, kinetics.reactions('VV'));
 Q_VV(i_out) = Q_VV_data;
end
if isKey(kinetics.reactions, 'Exch')
 if Exch_reactions(1).source == "Guerra95"
%                          N2 + O -> NO + N
  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, i1_N2)', ...
     Y(i_out, iO(1))', Y(i_out, kinetics.index{3})', ...
     Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reactions(1));
  Q_exch_N2_O(i_out) = Q_exch_NO_N_data;
%                          N2 + O <- NO + N
  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, i1_N2)', ...
     Y(i_out, iO(1))', Y(i_out, kinetics.index{3})', ...
     Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reactions(2));
  Q_exch_N_NO(i_out) = Q_exch_NO_N_data;
 else
%                          N2 + O -> NO + N
Exch_reaction_temp = Exch_reactions(1);
Exch_reaction_temp.reverse = false;
  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, iN2)', ...
     Y(i_out, iO)', Y(i_out, kinetics.index{3})', ...
     Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reaction_temp);
  Q_exch_N2_O(i_out) = Q_exch_NO_N_data;
%                          N2 + O <- NO + N
Exch_reaction_temp.direction_forward = false;
  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, iN2)', ...
     Y(i_out, iO)', Y(i_out, kinetics.index{3})', ...
     Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reaction_temp);
  Q_exch_N_NO(i_out) =  Q_exch_NO_N_data;

 end
end
end
c_p_N2 = c_p(N2, T);
c_p_O2 = c_p(O2, T);
c_p_total = 0.2*c_p_O2 + 0.8*c_p_N2;
Q_VT_Ks = Q_VT ./ (n_g/N_a) ./ c_p_total;
% Q_VT_wall_Ks = Q_VT_wall ./ (n_g/N_a) ./ c_p_total;
legend_str1(1) = "VT N_2-O, Pintassilgo2014";
legend_str2(1) = "VT N_2-O, code";
linewidth_Q = 2;
loglog(Pintassilgo2014_ag_Q_VT_N2_O(:, 1), ...
            Pintassilgo2014_ag_Q_VT_N2_O(:, 2), ...
                'color', [0 0 0.8], 'linewidth', linewidth_Q) %#ok<USENS>
hold on
if isKey(kinetics.reactions, 'Wall') && ...
                                    isKey(kinetics.reactions, 'Rec_wall')
 loglog(Pintassilgo2014_ag_Q_rec_O_wall(:, 1), ...
    Pintassilgo2014_ag_Q_rec_O_wall(:, 2), ...
                'color', [0.8 0 0], 'linewidth', linewidth_Q) %#ok<USENS>
 legend_str1 = [legend_str1, "O+wall rec, Pintassilgo2014"];
end
if isKey(kinetics.reactions, 'VV')
 loglog(Pintassilgo2014_ag_Q_VV_N2_N2(:, 1), ...
            Pintassilgo2014_ag_Q_VV_N2_N2(:, 2), ...
                'color', [0 0.7 0], 'linewidth', linewidth_Q) %#ok<USENS>
 legend_str1 = [legend_str1, "VV N2-N2, Pintassilgo2014"];
end
if isKey(kinetics.reactions, 'Exch')
 loglog(Pintassilgo2014_ag_Q_Zel_N2_O(:, 1), ...
            Pintassilgo2014_ag_Q_Zel_N2_O(:, 2), ...
               'color', [1 0.7 0], 'linewidth', linewidth_Q) %#ok<USENS>
 loglog(Pintassilgo2014_ag_Q_Zel_N_NO(:, 1), ...
            Pintassilgo2014_ag_Q_Zel_N_NO(:, 2), ...
               'color', [1 0.4 0], 'linewidth', linewidth_Q) %#ok<USENS>
 legend_str1 = [legend_str1, "Zel N2-O, Pintassilgo2014", ...
                                            "Zel N-NO, Pintassilgo2014"];
end
loglog(t_ag*1e3, Q_VT_Ks, ':', 'color', [0 0 0.8], ...
                                                'LineWidth', linewidth_Q)
if isKey(kinetics.reactions, 'Wall') && ...
                                    isKey(kinetics.reactions, 'Rec_wall')
 Q_rec_wall_Ks = Q_rec_wall ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_rec_wall_Ks, ':', 'color', [0.8 0 0], ...
                                                'LineWidth', linewidth_Q)
 legend_str2 = [legend_str2, "O+wall rec, code"];
end
if isKey(kinetics.reactions, 'VV')
 Q_VV_Ks = Q_VV ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_VV_Ks, ':', 'color', [0 0.7 0], ...
                                                'LineWidth', linewidth_Q)
 legend_str2 = [legend_str2, "VV N_2-N_2, code"];
end
if isKey(kinetics.reactions, 'Exch')
 Q_exch_N2_O_Ks = Q_exch_N2_O ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_exch_N2_O_Ks, ':', 'color', [1 0.7 0], ...
                                                'linewidth', linewidth_Q)
 Q_exch_N_NO_Ks = Q_exch_N_NO ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_exch_N_NO_Ks, ':', 'color', [1 0.4 0], ...
                                                'linewidth', linewidth_Q)
 legend_str2 = [legend_str2, "Zel N2-O, code", "Zel N-NO, code"];
end
%  loglog(t_ag*1e3, Q_VT_wall_Ks, 'linewidth', linewidth_Q)
legend([legend_str1 legend_str2], 'location', 'best')
title('Q_{in}')
xlim([6e-3 1e2])
ylim([6e0 1e5])
end

%% heating rates 2
if is_Q_plot2
 Q_data = zeros(length(Exch_reactions), length(t));
 for ind_exch = 1:length(Exch_reactions)
  reaction = Exch_reactions(ind_exch);
  IOM_M1 = IndexOfMolecules(reaction.particles(1));
  IOM_M2 = IndexOfMolecules(reaction.particles(2));
  IOM_M3 = IndexOfMolecules(reaction.particles(3));
  IOM_M4 = IndexOfMolecules(reaction.particles(4));
  M1 = kinetics.Ps{IOM_M1};
  M2 = kinetics.Ps{IOM_M2};
  M3 = kinetics.Ps{IOM_M3};
  M4 = kinetics.Ps{IOM_M4};
  indM1 = kinetics.index{IOM_M1};
  indM2 = kinetics.index{IOM_M2};
  indM3 = kinetics.index{IOM_M3};
  indM4 = kinetics.index{IOM_M4};
  for i_out = 1:length(t)
   switch length(reaction.particles)
    case 4
     [~, Q_temp] = R_exch({M1, M2, M3, M4}, Y(i_out, indM1)', ...
           Y(i_out, indM2)', Y(i_out, indM3)', Y(i_out, indM4)', ...
           T(i_out), reaction);
    case 5
     IOM_M5 = IndexOfMolecules(reaction.particles(5));
     indM5  = kinetics.index{IOM_M5};
     M5 = kinetics.Ps{IOM_M5};
     [~, Q_temp] = R_exch_23({M1, M2, M3, M4, M5}, ...
            Y(i_out, indM1)', Y(i_out, indM2)', Y(i_out, indM3)', ...
            Y(i_out, indM4)', Y(i_out, indM5)', T(i_out), reaction, 1);
    otherwise
     error("The number of particles is not 4 or 5.")
   end
   Q_data(ind_exch, i_out) = Q_temp;
  end
 end
end

Q_VT = zeros(1, length(t));
for i_out = 1:length(t)
 if isKey(kinetics.reactions, 'VT')
    [~, Q_VT_data] = R_VT(N2, Y(i_out, i1_N2)', O, ...
                Y(i_out, iO(1)), T(i_out), 1, kinetics.reactions('VT'));
    Q_VT(i_out) = Q_VT_data;
 end
end

Q_VT = Q_VT ./ (n_g'/N_a) ./ c_p_total';

Q_data = Q_data ./ (n_g'/N_a) ./ c_p_total';
figure
% hold on
title("Q_2")
l = {};
for ind_exch = 1:length(Exch_reactions)
    loglog(t_ag*1e3, Q_data(ind_exch, :), ...
                    'linewidth', 1.5, 'Color', Exch_colors(ind_exch, :));
    hold on
    l{ind_exch} = Exch_reactions(ind_exch).name; 
end
% xlim([0 20]);
legend(l);
xlim([6e-3 1e2])
ylim([6e0 1e5])
hold off
%% N2(A), N2(a'), N2(B) and N2(w) ag plot
if N2.num_elex_levels > 1
 i2_N2 = iN2(1 + N2.num_vibr_levels(1):...
                        N2.num_vibr_levels(1) + N2.num_vibr_levels(2));
 figure('Position', fsize) % N2(A), N2(a'), N2(B) and N2(w) ag plot
 loglog(Pintassilgo2014_ag_N2A(:, 1), Pintassilgo2014_ag_N2A(:, 2), ...
        'Color', [0.9 0 0], 'linewidth', 1.5, ...
                'DisplayName', "N2(A)/Ng, Pintassilgo2014") %#ok<USENS>
 hold on
 if N2.num_elex_levels > 2
  i3_N2 = iN2(1 + sum(N2.num_vibr_levels(1:2)):...
                        sum(N2.num_vibr_levels(1:3)));
  loglog(Pintassilgo2014_ag_N2B(:, 1), Pintassilgo2014_ag_N2B(:, 2), ...
            'color', [0 0.7 0], 'linewidth', 1.5, ...
                'DisplayName', 'N2(B)/Ng, Pintassilgo2014') %#ok<USENS>
 end
 loglog(t_ag*1e3, sum(Y(:, i2_N2), 2)./n_g, ':', ...
    'color', [0.9 0 0], 'linewidth', 1.5, 'DisplayName', "N2(A)/Ng, code")
 if N2.num_elex_levels > 2
  loglog(t_ag*1e3, sum(Y(:, i3_N2), 2)./n_g, ':', ...
    'color', [0 0.7 0], 'linewidth', 1.5, 'DisplayName', 'N2(B)/Ng, code')
 end
 % legend('N2(A)/Ng, Pintassilgo2014', ...'N2(B)/Ng, Pintassilgo2014', ...
 %     'N2(A)/Ng, code', ...'N2(B)/Ng, code', ...
 %                                                    'location', 'best')
 legend('location', 'best')
 xlabel('Afterglow time (ms)')
 xlim([1e-4 1.1e0])
 ylim([1e-8 1e-4])
 grid on
end

%% heating rates, K/s
% if N2.num_elex_levels > 2
%     figure('Position', fsize)
% iN2 = kinetics.index{1};
% iO2 = kinetics.index{2};
% % i1_N2 = iN2(1:N2.num_vibr_levels(1));
% % iO = kinetics.index{5};
% Q_N2B_O2 = zeros(length(t), 1);
% % Q_VT_wall = zeros(length(t), 1);
% % Q_rec_wall = zeros(length(t), 1);
% % Q_VV = zeros(length(t), 1);
% % Q_exch_N2_O = zeros(length(t), 1);
% % Q_exch_N_NO = zeros(length(t), 1);
% warning("recheck the number in Exch_reactions")
% for i_out = 1:length(t)
% [~, Q_N2B_O2_data] = R_exch_23({N2, O2, N2, O, O}, Y(i_out, iN2)', ...
%      Y(i_out, iO2)', Y(i_out, iN2)', Y(i_out, kinetics.index{5}(1))', ...
%      Y(i_out, kinetics.index{5}(1))', T(i_out), Exch_reactions(4));
% Q_N2B_O2(i_out) = Q_N2B_O2_data;
% % if isKey(kinetics.reactions, 'Wall') 
% % %  [~, Q_VT_wall_data] = ...
% % %                     R_VT_wall(N2, Y(i_out, i1_N2)', T(i_out), kinetics);
% % %  Q_VT_wall(i_out) = Q_VT_wall_data;
% %  if isKey(kinetics.reactions, 'Rec_wall')
% %   [~, Q_rec_wall_data] = R_rec_wall(O, Y(i_out, kinetics.index{5}), ...
% %                                                     T(i_out), kinetics);
% %   Q_rec_wall(i_out) = Q_rec_wall_data;
% %  end
% % end
% % if isKey(kinetics.reactions, 'VV')
% %  [~, Q_VV_data] = R_VV(N2, Y(i_out, i1_N2)', N2, Y(i_out, i1_N2)', ...
% %                              T(i_out), 1, 1, kinetics.reactions('VV'));
% %  Q_VV(i_out) = Q_VV_data;
% % end
% % if isKey(kinetics.reactions, 'Exch')
% % %                          N2 + O -> NO + N
% %  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, i1_N2)', ...
% %      Y(i_out, iO(1))', Y(i_out, kinetics.index{3})', ...
% %      Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reactions(1));
% %  Q_exch_N2_O(i_out) = Q_exch_NO_N_data;
% % %                          N2 + O <- NO + N
% %  [~, Q_exch_NO_N_data] = R_exch({N2, O, NO, N}, Y(i_out, i1_N2)', ...
% %      Y(i_out, iO(1))', Y(i_out, kinetics.index{3})', ...
% %      Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reactions(2));
% %  Q_exch_N_NO(i_out) = Q_exch_NO_N_data;
% % end
% end
% Q_N2B_O2_Ks = Q_N2B_O2 ./ (n_g/N_a) ./ c_p_total;
% clear legend_str1 legend_str2
% legend_str1(1) = "N_2(B)+O_2, Pintassilgo2014";
% legend_str2(1) = "N_2(B)+O_2, code";
% loglog(Pintassilgo2014_ag_Q_N2B_O2(:, 1), ...
%             Pintassilgo2014_ag_Q_N2B_O2(:, 2), ...
%                         'color', [0 0 0.8], 'linewidth', 1.5) %#ok<USENS>
% hold on
% % if isKey(kinetics.reactions, 'Wall') && ...
% %                                     isKey(kinetics.reactions, 'Rec_wall')
% %  loglog(Pintassilgo2014_ag_Q_rec_O_wall(:, 1), ...
% %     Pintassilgo2014_ag_Q_rec_O_wall(:, 2), ...
% %                         'color', [0.8 0 0], 'linewidth', 1.5) %#ok<USENS>
% %  legend_str1 = [legend_str1, "O+wall rec, Pintassilgo2014"];
% % end
% % if isKey(kinetics.reactions, 'VV')
% %  loglog(Pintassilgo2014_ag_Q_VV_N2_N2(:, 1), ...
% %             Pintassilgo2014_ag_Q_VV_N2_N2(:, 2), ...
% %                         'color', [0 0.7 0], 'linewidth', 1.5) %#ok<USENS>
% %  legend_str1 = [legend_str1, "VV N2-N2, Pintassilgo2014"];
% % end
% % if isKey(kinetics.reactions, 'Exch')
% %  loglog(Pintassilgo2014_ag_Q_Zel_N2_O(:, 1), ...
% %             Pintassilgo2014_ag_Q_Zel_N2_O(:, 2), ...
% %                        'color', [1 0.7 0], 'linewidth', 1.5) %#ok<USENS>
% %  loglog(Pintassilgo2014_ag_Q_Zel_N_NO(:, 1), ...
% %             Pintassilgo2014_ag_Q_Zel_N_NO(:, 2), ...
% %                        'color', [1 0.4 0], 'linewidth', 1.5) %#ok<USENS>
% %  legend_str1 = [legend_str1, "Zel N2-O, Pintassilgo2014", ...
% %                                             "Zel N-NO, Pintassilgo2014"];
% % end
% loglog(t_ag*1e3, Q_N2B_O2_Ks, ':', 'color', [0 0 0.8], 'linewidth', 1.5)
% % if isKey(kinetics.reactions, 'Wall') && ...
% %                                     isKey(kinetics.reactions, 'Rec_wall')
% %  Q_rec_wall_Ks = Q_rec_wall ./ (n_g/N_a) ./ c_p_total;
% %  loglog(t_ag*1e3, Q_rec_wall_Ks, ':', 'color', [0.8 0 0], 'linewidth',1.5)
% %  legend_str2 = [legend_str2, "O+wall rec, code"];
% % end
% % if isKey(kinetics.reactions, 'VV')
% %  Q_VV_Ks = Q_VV ./ (n_g/N_a) ./ c_p_total;
% %  loglog(t_ag*1e3, Q_VV_Ks, ':', 'color', [0 0.7 0], 'linewidth', 1.5)
% %  legend_str2 = [legend_str2, "VV N_2-N_2, code"];
% % end
% % if isKey(kinetics.reactions, 'Exch')
% %  Q_exch_N2_O_Ks = Q_exch_N2_O ./ (n_g/N_a) ./ c_p_total;
% %  loglog(t_ag*1e3, Q_exch_N2_O_Ks, ':', 'color', [1 0.7 0], ...
% %                                                         'linewidth', 1.5)
% %  Q_exch_N_NO_Ks = Q_exch_N_NO ./ (n_g/N_a) ./ c_p_total;
% %  loglog(t_ag*1e3, Q_exch_N_NO_Ks, ':', 'color', [1 0.4 0], ...
% %                                                         'linewidth', 1.5)
% %  legend_str2 = [legend_str2, "Zel N2-O, code", "Zel N-NO, code"];
% % end
% % %  loglog(t_ag*1e3, Q_VT_wall_Ks, 'linewidth', 1.5)
% legend([legend_str1 legend_str2], 'location', 'best')
% title('Q_{in} with el. excitation')
% xlim([6e-3 1e2])
% ylim([6e0 1e5])
% end

%% N2+ ag plot
if length(kinetics.Ps) > 5
 figure
 loglog(t_ag*1e3, Y(:, kinetics.index{IndexOfMolecules("N2+")})./n_g, ...
      ':', 'color', [0.9 0 0], 'linewidth', 1.5, 'DisplayName', 'N2+/n_g')
 if length(kinetics.Ps) > 6
  hold on
  loglog(t_ag*1e3, Y(:, kinetics.index{IndexOfMolecules("O2+")})./n_g, ...
      ':', 'color', [0 0.7 0], 'linewidth', 1.5, 'DisplayName', 'O2+/n_g')
  if isKey(kinetics.reactions, 'free_e')
   loglog(t_ag*1e3, Y(:, end-2)./n_g, ...
         ':', 'color', [0 0 0], 'linewidth', 1.5, 'DisplayName', 'e-/n_g')
  end
 end
 legend('location', 'best')
 xlabel('Afterglow time (ms)')
 xlim([1e-4 1.1e0])
 ylim([1e-8 1e-4])
end

rmpath('../src/')
end