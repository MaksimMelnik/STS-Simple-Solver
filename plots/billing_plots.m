load('../data/particles.mat', 'N2', 'O2', 'N', 'O', 'NO');
addpath('../src/')

for i = 4%[1 2] % [1 2 3 4]   choosing reactions
                        % 1 is for N2(i) + N2 -> N2(i-1) + N2
                        % 2 is for N2(i) + N  -> N2(i-1) + N
                        % 3 is for O2(i) + O2 -> O2(i-1) + O2
                        % 4 is for O2(i) + O  -> O2(i-1) + O
    if i == 1
        figure
        semilogy(kvt_ssh(1000, N2, N2, 1, 1));
        hold on
        semilogy(kvt_fho_old(1000, N2, N2, 1));
        hold on
        semilogy(kvt_billing(1000, N2, N2, 1, 1));
        hold off
        legend('N2+N2 SSH', 'N2+N2 FHO_{OLD}', 'N2+N2 BILLING');
        title('Comparing of VT-excharge rate coefficients in nitrogen')
    elseif i == 2
        figure
        semilogy(kvt_ssh(1000, N2, N, 1, 1));
        hold on
        semilogy(kvt_fho_old(1000, N2, N, 1));
        hold on
        semilogy(kvt_billing(1000, N2, N, 1, 1));
        hold off
        legend('N2+N SSH', 'N2+N FHO_OLD', 'N2+N BILLING');
        title('Comparing of VT-excharge rate coefficients in nitrogen')
    elseif i == 3
        figure
        semilogy(kvt_ssh(1000, O2, O2, 1, 1));
        hold on
        semilogy(kvt_fho_old(1000, O2, O2, 1));
        hold on
        semilogy(kvt_billing(1000, O2, O2, 1, 1));
        hold off
        legend('O2+O2 SSH', 'O2+O2 FHO_OLD', 'O2+O2 BILLING');
        title('Comparing of VT-excharge rate coefficients in oxygen')
    elseif i == 4
        figure
        semilogy(kvt_ssh(1000, O2, O, 1, 1));
        hold on
        semilogy(kvt_fho_old(1000, O2, O, 1));
        hold on
        semilogy(kvt_billing(1000, O2, O, 1, 1));
        hold off
        legend('O2+O SSH', 'O2+O FHO_OLD', 'O2+O BILLING');
        title('Comparing of VT-excharge rate coefficients in oxygen')
    elseif i == 5
        T = 6000;
        figure
        semilogy(kvt_billing(T, N2, N2, 1, 1));
        hold on
        semilogy(kvt_billing(T, N2, N, 1, 1));
        legend('N2+N2 B', 'N2+N B', 'location', 'best');
        title('Comparing of VT-excharge rate coefficients in nitrogen')
        ylim([5e-19 3e-15])

    end
end