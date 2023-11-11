load('../STS-Simple-Solver/data/particles.mat', 'N2', 'O2', 'N', 'O', 'NO');

for i = [1 2] % [1 2 3 4]   choosing reactions
                        % 1 is for N2(i) + N2 -> N2(i-1) + N2
                        % 2 is for N2(i) + N  -> N2(i-1) + N
                        % 3 is for O2(i) + O2 -> O2(i-1) + O2
                        % 4 is for O2(i) + O  -> O2(i-1) + O
    if i == 1
        figure
        plot(kvt_ssh(1000, N2, N2, 1, 1));
        hold on
        plot(kvt_fho_old(1000, N2, N2, 1));
        hold on
        plot(kvt_billing(1000, N2, N2, 1, 1));
        hold off
        legend('N2+N2 SSH', 'N2+N2 FHO_OLD', 'N2+N2 BILLING');
        title('Comparing of VT-excharge rate coefficients in nitrogen')
    elseif i == 2
        figure
        plot(kvt_ssh(1000, N2, N, 1, 1));
        hold on
        plot(kvt_fho_old(1000, N2, N, 1));
        hold on
        plot(kvt_billing(1000, N2, N, 1, 1));
        hold off
        legend('N2+N SSH', 'N2+N FHO_OLD', 'N2+N BILLING');
        title('Comparing of VT-excharge rate coefficients in nitrogen')
    elseif i == 3
        figure
        plot(kvt_ssh(1000, O2, O2, 1, 1));
        hold on
        plot(kvt_fho_old(1000, O2, O2, 1));
        hold on
        plot(kvt_billing(1000, O2, O2, 1, 1));
        hold off
        legend('O2+O2 SSH', 'O2+O2 FHO_OLD', 'O2+O2 BILLING');
        title('Comparing of VT-excharge rate coefficients in oxygen')
    elseif i == 4
        figure
        plot(kvt_ssh(1000, O2, O, 1, 1));
        hold on
        plot(kvt_fho_old(1000, O2, O, 1));
        hold on
        plot(kvt_billing(1000, O2, O, 1, 1));
        hold off
        legend('O2+O SSH', 'O2+O FHO_OLD', 'O2+O BILLING');
        title('Comparing of VT-excharge rate coefficients in oxygen')

    end
end