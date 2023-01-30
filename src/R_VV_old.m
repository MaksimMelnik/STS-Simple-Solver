function out=R_VV_old(M1, n_m, M2, n_M2, T, ind_e1, ind_e2, model)
% Universal function for relaxation parts R_VV calculation.
% M1 is the first molecule, n_m is the array of M1's number
% densities (n_i), M2 is the collision partner, n_M2 is the number density
% of M2, T is the gas temperature, ind_e1 is the electronic level of M1,
% ind_e2 is the electronic level of M2, model (optional) is the VT model.
% 20.01.2023 Maksim Melnik

if nargin<8
    model='SSH';
end

k = 1.380649e-23;

switch model
    case 'FHO'
        error('FHO VV is not implemented yet')
%         k_down=kvv_fho_old(T, M1, M2, ind_e1)';
    case 'SSH'
        k_down=kvv_ssh(T, M1, M2, ind_e1, ind_e2, 1);
end
k_up=k_down.*exp((M1.ev_i{ind_e1}(1:end-1)'+M2.ev_i{ind_e2}(2:end)...
                 -M1.ev_i{ind_e1}(2:end)'-M2.ev_i{ind_e2}(1:end-1))/k/T);
core1=2:M1.num_vibr_levels(ind_e1)-1; % c 1 по n-1 уровень
core2=2:M2.num_vibr_levels(ind_e2)-1; % c 1 по n-1 уровень
R=zeros(M1.num_vibr_levels(ind_e1), M2.num_vibr_levels(ind_e2));
    % заполняем серединку
R(core1, core2)=k_down(core1, core2-1).* n_m(core1+1) .* n_M2(core2-1)'...
                 + k_up(core1-1, core2).*n_m(core1-1).* n_M2(core2+1)' ...
                   - k_down(core1-1, core2).*n_m(core1).* n_M2(core2)' ...
                        - k_up(core1, core2-1).*n_m(core1).*n_M2(core2)';
    % заполняем для первого и последнего уровней
R(1, 2:end)=k_down(1, :).* n_m(2) .* n_M2(1:end-1)' ...
                                - k_up(1, :) .* n_m(1) .* n_M2(2:end)';
R(end, 1:end-1)=k_up(end, :).*n_m(end-1).* n_M2(2:end)' ...
                            - k_down(end, :).*n_m(end).* n_M2(1:end-1)';
R(2:end, 1)= k_up(:, 1).*n_m(1:end-1).* n_M2(2) ...
                                    - k_down(:, 1).*n_m(2:end).* n_M2(1);
R(1:end-1, end)=k_down(:, end).* n_m(2:end) .* n_M2(end-1) ...
                                - k_up(:, end).*n_m(1:end-1).*n_M2(end);
out=R;
end