function out=R_VT_old(M, n_m, P, n_p, T, ind_e)
% (T, n1, n3, n4, n5)
% универсальная функция расчёта релаксационных членов R_VT, но пока только
% для CO, потому что там k_VT для CO только подходит.
% 11.10.2020
k = 1.380649e-23; 
core=2:M.num_vibr_levels(ind_e)-1; %c 1 по n-1 уровень

R=zeros(M.num_vibr_levels(ind_e), 1);
switch P.name
    case 'CO'
        ind_p=1;
    case 'C'
        ind_p=2;
    case 'O'
        ind_p=3;
    case 'Ar'
        ind_p=4;
    case 'C2'
        ind_p=5;
end
num_v_l=M.num_vibr_levels(ind_e);
k_down=kvt_fho_old(T, M, P, ind_e)';
k_up=k_down.*exp((M.ev_i{ind_e}(1:end-1)-M.ev_i{ind_e}(2:end))/k/T)';
% disp(Kvt(1))
    % заполняем серединку
R(core)=n_p*(k_down(core).*n_m(core+1)+k_up(core-1).*n_m(core-1)...
                       -k_down(core-1).*n_m(core)-k_up(core).*n_m(core));
    % заполняем для первого и последнего уровней
R(1)=n_p*(k_down(1)*n_m(2)-k_up(1)*n_m(1));
R(M.num_vibr_levels(ind_e))=...
  n_p*(k_up(M.num_vibr_levels(ind_e)-1)*n_m(M.num_vibr_levels(ind_e)-1)...
    -k_down(M.num_vibr_levels(ind_e)-1)*n_m(M.num_vibr_levels(ind_e)));
out=R;
end