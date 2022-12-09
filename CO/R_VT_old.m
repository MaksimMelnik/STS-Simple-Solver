function out=R_VT_old(M1, n_m, M2, n_M2, T, ind_e, model)
% (T, n1, n3, n4, n5)
% универсальная функция расчёта релаксационных членов R_VT, но пока только
% для CO, потому что там k_VT для CO только подходит.
% 11.10.2020
% M1 is the first molecule, n_m is the array of M1's number
% densities (n_i), M2 is the collision partner, n_M2 is the number density
% of M2, T is the gas temperature, ind_e is the electronic level of M1,
% model (optional) is the VT model.

if nargin<7
    model='FHO';
end

k = 1.380649e-23;

switch model
    case 'FHO'
        k_down=kvt_fho_old(T, M1, M2, ind_e)';
    case 'SSH'
        k_down=kvt_ssh(T, M1, M2, ind_e, 1)';
end
core=2:M1.num_vibr_levels(ind_e)-1; % c 1 по n-1 уровень
R=zeros(M1.num_vibr_levels(ind_e), 1);
k_up=k_down.*exp((M1.ev_i{ind_e}(1:end-1)-M1.ev_i{ind_e}(2:end))/k/T)';
    % заполняем серединку
R(core)=n_M2*(k_down(core).*n_m(core+1)+k_up(core-1).*n_m(core-1)...
                       -k_down(core-1).*n_m(core)-k_up(core).*n_m(core));
    % заполняем для первого и последнего уровней
R(1)=n_M2*(k_down(1)*n_m(2)-k_up(1)*n_m(1));
R(M1.num_vibr_levels(ind_e))=...
  n_M2*(k_up(M1.num_vibr_levels(ind_e)-1)*n_m(M1.num_vibr_levels(ind_e)-1)...
    -k_down(M1.num_vibr_levels(ind_e)-1)*n_m(M1.num_vibr_levels(ind_e)));
out=R;
end