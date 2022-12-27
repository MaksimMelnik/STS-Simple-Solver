function out=R_VE_m(M, n_m, n_ex, n_p, T, ind_e)
% функция расчёта R_VE для одного электронного состояния, прыжки со всех
% колебательных уровней основного на все колебательные возбуждённого (в
% рамках резонанса, конечно). 
% Модификация, основанная на коде без колебательных уровней возбуждённого.
% Only for CO.
% 29.07.2021
k=1.380649e-23; h = 6.626070041e-34; c=299792458;
k_VE_up=k_VE_m(M, T, ind_e);

Theta_r_M = M.Be(1)*h*c/k;
Zrot_M = T./(M.sigma.*Theta_r_M);
Theta_r_ex= M.Be(ind_e)*h*c/k;
Zrot_ex= T./(M.sigma.*Theta_r_ex);
dE=M.e_E(ind_e)+M.ev_0(ind_e)+M.ev_i{ind_e}-(M.ev_i{1}+M.ev_0(1))';
k_VE_down=k_VE_up*M.s_e(1)/M.s_e(ind_e).*Zrot_M/Zrot_ex.*exp(dE/k/T);
R=n_p*(n_ex'.*k_VE_down-n_m.*k_VE_up);
out=R;
end