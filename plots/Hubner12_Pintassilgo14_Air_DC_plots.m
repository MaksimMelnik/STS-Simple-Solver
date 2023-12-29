function Hubner12_Pintassilgo14_Air_DC_plots(data)
% Comparison plots for simulations of Hubner experiment [1] with previous
% simulations by Pintassilgo [2].
% [1] M Hubner et al Meas. Sci. Technol. 23 (2012) 115602.
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.
% 26.07.2023 Maksim Melnik

is_T_Tv_plot = true;
is_Q_plot = true;
is_n_NO_N_O = true;
is_VDF = true;

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
if isKey(kinetics.reactions, 'Exch')
    Exch_reactions = kinetics.reactions("Exch");
    IndexOfMolecules = kinetics.IndexOfMolecules;
end

n_g=sum(Y(:, 1:end-1), 2);
n_N2=sum(Y(:, kinetics.index{1}), 2);
t_ag=t-0.005;
fsize = [200 50 900 550];

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
plot(Hubner_2012_T(:, 1), Hubner_2012_T(:, 2), 'sq', t*1e3, T, ...
                        t*1e3, Tv, '-.', 'linewidth', 1.5) %#ok<USENS>
% errorbar T Hubner +- 40 K
legend('T_{exp}, Hubner 2012', 'T', 'Tv', 'location', 'best')
xlabel('t, ms')
xlim([-2 14])
ylim([250 620])

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
    figure('Position', fsize)
time_ind0=1;
[~, time_ind1]   = min( abs(t_ag*1e3 - 1) );
[~, time_ind10]  = min( abs(t_ag*1e3 - 10) );
% time_ind100=length(Y(:, 1));
[~, time_ind100] = min( abs(t_ag*1e3 - 100) );
lvls4plot=0:length(Y(1, kinetics.index{1}))-1;
semilogy(Pintassilgo2014_N2_VDF_post_DC(:,1), ...
            Pintassilgo2014_N2_VDF_post_DC(:,2), ...
                           'color', [1 0 0] ,'linewidth', 1.5) %#ok<USENS>
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
[~, Q_VT_data] = R_VT(N2, Y(i_out, i1_N2)', O, ...
                Y(i_out, iO(1)), T(i_out), 1, kinetics.reactions('VT'));
Q_VT(i_out) = Q_VT_data;
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
loglog(Pintassilgo2014_ag_Q_VT_N2_O(:, 1), ...
            Pintassilgo2014_ag_Q_VT_N2_O(:, 2), ...
                        'color', [0 0 0.8], 'linewidth', 1.5) %#ok<USENS>
hold on
if isKey(kinetics.reactions, 'Wall') && ...
                                    isKey(kinetics.reactions, 'Rec_wall')
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
 loglog(Pintassilgo2014_ag_Q_Zel_N_NO(:, 1), ...
            Pintassilgo2014_ag_Q_Zel_N_NO(:, 2), ...
                       'color', [1 0.4 0], 'linewidth', 1.5) %#ok<USENS>
 legend_str1 = [legend_str1, "Zel N2-O, Pintassilgo2014", ...
                                            "Zel N-NO, Pintassilgo2014"];
end
loglog(t_ag*1e3, Q_VT_Ks, ':', 'color', [0 0 0.8], 'linewidth', 1.5)
if isKey(kinetics.reactions, 'Wall') && ...
                                    isKey(kinetics.reactions, 'Rec_wall')
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
 Q_exch_N_NO_Ks = Q_exch_N_NO ./ (n_g/N_a) ./ c_p_total;
 loglog(t_ag*1e3, Q_exch_N_NO_Ks, ':', 'color', [1 0.4 0], ...
                                                        'linewidth', 1.5)
 legend_str2 = [legend_str2, "Zel N2-O, code", "Zel N-NO, code"];
end
%  loglog(t_ag*1e3, Q_VT_wall_Ks, 'linewidth', 1.5)
legend([legend_str1 legend_str2], 'location', 'best')
title('Q_{in}')
xlim([6e-3 1e2])
ylim([6e0 1e5])
end
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
                                ':', 'color', [0.9 0 0], 'linewidth', 1.5)
 if length(kinetics.Ps) > 6
  hold on
  loglog(t_ag*1e3, Y(:, kinetics.index{IndexOfMolecules("O2+")})./n_g, ...
                                ':', 'color', [0 0.7 0], 'linewidth', 1.5)
  % loglog(t_ag*1e3, Y(:, end-1)./n_g, ...
                                % ':', 'color', [0 0 0], 'linewidth', 1.5)
 end
 legend('N2+/n_g', 'O2+/n_g', 'e-/n_g', 'location', 'best')
 xlabel('Afterglow time (ms)')
 xlim([1e-4 1.1e0])
 ylim([1e-8 1e-4])
end

rmpath('../src/')
end