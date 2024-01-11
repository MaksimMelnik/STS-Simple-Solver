function out=density_f_exc(T, n1, M)
% Равновесное распределение CO по колебательным уровням с учётом 
% возбуждённых состояний.
% 29.07.2021
k = 1.380649e-23;
if M.fr_deg_c > 3
 if M.num_elex_levels==1
  Zv=sum(exp(-M.ev_i{1}/T/k));
  n=n1/Zv.*exp(-M.ev_i{1}/T/k)';
 else
    if M.num_elex_levels ~= 3
        error("universal distribution is still not implemented")
    end
  Zel=M.s_e*exp(-M.e_E/k/T)';
  Zvibr=[ sum(exp(-(M.ev_i{1}+M.ev_0(1))/k/T));
         sum(exp(-(M.ev_i{2}+M.ev_0(2))/k/T));
         sum(exp(-(M.ev_i{3}+M.ev_0(3))/k/T))];


  n=n1/Zel*[  M.s_e(1)*exp(-(M.ev_i{1}+M.ev_0(1)+M.e_E(1))/k/T)'/Zvibr(1);
             M.s_e(2)*exp(-(M.ev_i{2}+M.ev_0(2)+M.e_E(2))/k/T)'/Zvibr(2);
             M.s_e(3)*exp(-(M.ev_i{3}+M.ev_0(3)+M.e_E(3))/k/T)'/Zvibr(3)];
 end
else
    if M.num_elex_levels > 1
        error("Equilibrium distribution for atoms is not implemented yet")
    end
    n = n1;
end
out=n;
end