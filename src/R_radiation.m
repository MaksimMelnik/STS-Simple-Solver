function [R, Q] = R_radiation(P1, n1, P2, n2)
% The test function for calculation of relaxation term R for radiation.
% currently only for N2(B) → N2(A) + hv
% source: Pintassilgo2009
% 25.10.2024 Maksim Melnik

h = 6.62607015e-34;     % Plank constant, J*sec
v_radiation = 2e5;      % s-1, from Pintassilgo2009
    % stupid simple model
if sum(n1) > 0
    n1_d = n1 / sum(n1);
else
    n1_d = n1 * 0;
end
if sum(n2) > 0
    n2_d = n2 / sum(n2);
else
    n2_d = n2 * 0;
end
v_matrix = n1_d .* n2_d';
R = - n1 .* v_matrix * v_radiation;
ind_e1 = 3;
ind_e2 = 2;
dE = h*v_radiation + P2.e_E(ind_e2) + P2.ev_0(ind_e2) + P2.ev_i{ind_e2}...
                    - P1.e_E(ind_e1) - P1.ev_0(ind_e1) - P1.ev_i{ind_e1}';
Q = sum(R .* dE, 'all');
end