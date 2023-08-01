function Hubner12_Pintassilgo14_Air_DC_plots(data)
% Comparison plots for simulations of Hubner experiment [1] with previous
% simulations by Pintassilgo [2].
% [1] M Hubner et al Meas. Sci. Technol. 23 (2012) 115602.
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.
% 26.07.2023 Maksim Melnik

    % constants
N_a = 6.02214076e23;              % Avogadro constant
    % loading data from papers
load('../data/for comparison/Hubner2012_and_Pintassilgo2014.mat' ...
                                                            ) %#ok<LOAD>
addpath('../src/')

    % taking required data
Y  = data.Y;
t  = data.t;
T  = data.T;
Tv = data.Tv;
kinetics = data.kinetics;
N2 = kinetics.Ps{1};
O2 = kinetics.Ps{2};
NO = kinetics.Ps{3};
N  = kinetics.Ps{4};
O  = kinetics.Ps{5};
Exch_reactions = kinetics.reactions("Exch");

n_g=sum(Y(:, 1:end-1), 2);
n_N2=sum(Y(:, kinetics.index{1}), 2);
t_ag=t-0.005;
fsize = [200 50 900 550];


    figure  % T and Tv plot
plot(t*1e3, T, t*1e3, Tv, '-.', 'linewidth', 1.5)
legend('T', 'Tv', 'location', 'best')
xlabel('t, ms')
% xlim([6e-2 1.5e-1])


    figure  % T vs exp plot
plot(Hubner_2012_T(:, 1), Hubner_2012_T(:, 2), 'sq', t*1e3, T, ...
                        t*1e3, Tv, '-.', 'linewidth', 1.5) %#ok<USENS>
% errorbar T Hubner +- 40 K
legend('T_{exp}, Hubner 2012', 'T', 'Tv', 'location', 'best')
xlabel('t, ms')
xlim([-2 14])
ylim([250 620])


    figure('Position', fsize)  % N, O and NO ag plot
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


    figure('Position', fsize) % VDF ag plot
time_ind0=1;
[~, time_ind1] = min( abs(t_ag*1e3 - 1) );
[~, time_ind10] = min( abs(t_ag*1e3 - 10) );
time_ind100=length(Y(:, 1));
lvls4plot=0:length(Y(1, kinetics.index{1}))-1;
semilogy(Pintassilgo2014_N2_VDF_post_DC(:,1), ...
            Pintassilgo2014_N2_VDF_post_DC(:,2), ...
                                'color', [1 0 0] ,'linewidth', 1.5)
hold on
semilogy(Pintassilgo2014_N2_VDF_ag_1ms(:,1), ...
                            Pintassilgo2014_N2_VDF_ag_1ms(:,2), ...
                'color', [0.75 0.2 0.25] ,'linewidth', 1.5) %#ok<USENS>
semilogy(Pintassilgo2014_N2_VDF_ag_10ms(:,1), ...
                            Pintassilgo2014_N2_VDF_ag_10ms(:,2), ...
                'color', [0.25 0.4 0.75] ,'linewidth', 1.5) %#ok<USENS>
semilogy(Pintassilgo2014_N2_VDF_ag_100ms(:,1), ...
                        Pintassilgo2014_N2_VDF_ag_100ms(:,2), ...
                        'color', [0 0.6 1] ,'linewidth', 1.5) %#ok<USENS>
semilogy(lvls4plot, Y(time_ind0, kinetics.index{1})/n_N2(time_ind0), ...
                            '--', 'color', [1 0 0] ,'linewidth', 2.5)
semilogy(lvls4plot, Y(time_ind1, kinetics.index{1})/n_N2(time_ind1), ...
                        '--', 'color', [0.75 0.2 0.25] ,'linewidth', 2.5)
semilogy(lvls4plot, Y(time_ind10, kinetics.index{1})/n_N2(time_ind10), ...
                        '--', 'color', [0.25 0.4 0.75] ,'linewidth', 2.5)
semilogy(lvls4plot, ...
            Y(time_ind100, kinetics.index{1})/n_N2(time_ind100), ...
                            '--', 'color', [0 0.6 1] ,'linewidth', 2.5)
 legend('0.1 ms, Pintassilgo2014', '1 ms, Pintassilgo2014', ...
            '10 ms, Pintassilgo2014', '100 ms, Pintassilgo2014', ...
            '0 ms, Maksim', '1 ms, Maksim', '10 ms, Maksim', ...
                                    '100 ms, Maksim', 'location', 'best')
xlim([0 30])
ylim([1e-6 1])


    figure('Position', fsize)  % heating rates, K/s
iN2 = kinetics.index{1};
i1_N2 = iN2(1:N2.num_vibr_levels(1));
iO = kinetics.index{5};
Q_VT = zeros(length(t), 1);
Q_VT_wall = zeros(length(t), 1);
Q_rec_wall = zeros(length(t), 1);
Q_VV = zeros(length(t), 1);
Q_exch_N2_O = zeros(length(t), 1);
for i_out = 1:length(t)
[~, Q_VT_data] = R_VT(N2, Y(i_out, i1_N2)', O, ...
                Y(i_out, iO(1)), T(i_out), 1, kinetics.reactions('VT'));
Q_VT(i_out) = Q_VT_data;
if isKey(kinetics.reactions, 'Wall') 
%  [~, Q_VT_wall_data] = ...
%                     R_VT_wall(N2, Y(i_out, i1_N2)', T(i_out), kinetics);
%  Q_VT_wall(i_out) = Q_VT_wall_data;
 if isKey(kinetics.reactions, 'Diss')
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
%                          N2 + O -> NO + N
 [~, Q_exch_NO_N_data] = R_exch(N2, O, NO, N, Y(i_out, i1_N2)', ...
     Y(i_out, iO(1))', Y(i_out, kinetics.index{3})', ...
     Y(i_out, kinetics.index{4}(1))', T(i_out), Exch_reactions(1));
 Q_exch_N2_O(i_out) = Q_exch_NO_N_data;
end
end
c_p_N2 = c_p(N2, T);
c_p_O2 = c_p(O2, T);
c_p_total = 0.2*c_p_O2 + 0.8*c_p_N2;
Q_VT_Ks = Q_VT ./ (n_g/N_a) ./ c_p_total;
% Q_VT_wall_Ks = Q_VT_wall ./ (n_g/N_a) ./ c_p_total;
legend_str1 = ["VT N_2-O, Pintassilgo2014"];
legend_str2 = ["VT N_2-O, code"];
loglog(Pintassilgo2014_ag_Q_VT_N2_O(:, 1), ...
            Pintassilgo2014_ag_Q_VT_N2_O(:, 2), ...
                        'color', [0 0 0.8], 'linewidth', 1.5) %#ok<USENS>
hold on
if isKey(kinetics.reactions, 'Wall') && isKey(kinetics.reactions, 'Diss')
 loglog(Pintassilgo2014_ag_Q_rec_O_wall(:, 1), ...
    Pintassilgo2014_ag_Q_rec_O_wall(:, 2), ...
                        'color', [0.8 0 0], 'linewidth', 1.5) %#ok<USENS>
 legend_str1 = [legend_str1, "O+wall rec, Pintassilgo2014"];
end
if isKey(kinetics.reactions, 'VV')
 loglog(Pintassilgo2014_ag_Q_VV_N2_N2(:, 1), ...
            Pintassilgo2014_ag_Q_VV_N2_N2(:, 2), ...
                        'color', [0 0.7 0], 'linewidth', 1.5) %#ok<USENS>
 legend_str1 = [legend_str1, "VV N2-N2, Pintassilgo2014"];
end
if isKey(kinetics.reactions, 'Exch')
 loglog(Pintassilgo2014_ag_Q_Zel_N2_O(:, 1), ...
            Pintassilgo2014_ag_Q_Zel_N2_O(:, 2), ...
                       'color', [1 0.7 0], 'linewidth', 1.5) %#ok<USENS>
 legend_str1 = [legend_str1, "Zel N2-O, Pintassilgo2014"];
end
loglog(t_ag*1e3, Q_VT_Ks, ':', 'color', [0 0 0.8], 'linewidth', 1.5)
if isKey(kinetics.reactions, 'Wall') && isKey(kinetics.reactions, 'Diss')
 Q_rec_wall_Ks = Q_rec_wall ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_rec_wall_Ks, ':', 'color', [0.8 0 0], 'linewidth',1.5)
 legend_str2 = [legend_str2, "O+wall rec, code"];
end
if isKey(kinetics.reactions, 'VV')
 Q_VV_Ks = Q_VV ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_VV_Ks, ':', 'color', [0 0.7 0], 'linewidth', 1.5)
 legend_str2 = [legend_str2, "VV N_2-N_2, code"];
end
if isKey(kinetics.reactions, 'Exch')
 Q_exch_N2_O_Ks = Q_exch_N2_O ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_exch_N2_O_Ks, ':', 'color', [1 0.7 0], ...
                                                        'linewidth', 1.5) 
 legend_str2 = [legend_str2, "Zel NO-N, code"];
end
%  loglog(t_ag*1e3, Q_VT_wall_Ks, 'linewidth', 1.5)
legend([legend_str1 legend_str2], 'location', 'best')
title('Q_{in}')
xlim([6e-3 1e2])
ylim([6e0 1e5])

rmpath('../src/')
end