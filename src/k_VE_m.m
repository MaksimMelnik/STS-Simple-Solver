function out=k_VE_m(M, T, ind_e)
% расчёт коэффициентов VE-обменов CO(X)->CO(a3П, A1П) в соответствие с
% данными в дипломе Мишиной.
% fix в матричной форме для произвольных i и j
% 13.07.2021
k=1.380649e-23; h = 6.626070041e-34; c=299792458;
dE=M.e_E(ind_e)+M.ev_0(ind_e)+M.ev_i{ind_e} - (M.ev_i{1}+M.ev_0(1))';
we=max(M.we(1), M.we(ind_e));
out=M.S_VE(ind_e)*...
              exp(-(dE/(h*c)/(M.C_VE(ind_e)*we)).^2) .*exp(-dE/(2*k*T));
end