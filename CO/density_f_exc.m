function out=density_f_exc(T, n1, M)
% Равновесное распределение CO по колебательным уровням с учётом 
% возбуждённых состояний.
% 29.07.2021
k = 1.380649e-23; c = 299792458; h = 6.626070041e-34;
% n=zeros(sum(M.num_vibr_levels),1);
% Theta_r = M.Be*h*c/k;
% Z_rot = T./(M.sigma.*Theta_r);
% Zint=0;
    % пробегаем по всем электронным состояниям
% for i=1:M.num_elex_levels
%     
%  Zint=Zint+M.s_e(i)...
%     *sum(exp(-(M.e_E(i)+M.e_0(i)+M.e_i(i,1:M.num_vibr_levels(i)))/T/k));
%  n(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
%     n1*M.s_e(i)...*Z_rot(i)...
%         *exp(-(M.e_E(i)+M.e_0(i)+M.e_i(i,1:M.num_vibr_levels(i)))/T/k)';
%  n(1+sum(M.num_vibr_levels(1:i-1)):sum(M.num_vibr_levels(1:i)))=...
%     n1*M.s_e(i)...*Z_rot(i)...
%         *exp(-(M.e_E(i)+M.e_0(i)+M.e_i(i,1:M.num_vibr_levels(i)))/T/k)';
% end
% out=n/Zint;

Zel=M.s_e*exp(-M.e_E/k/T)';
Zvibr=[ sum(exp(-(M.ev_i{1}+M.ev_0(1))/k/T));
        sum(exp(-(M.ev_i{2}+M.ev_0(2))/k/T));
        sum(exp(-(M.ev_i{3}+M.ev_0(3))/k/T))];


n=n1/Zel*[  M.s_e(1)*exp(-(M.ev_i{1}+M.ev_0(1)+M.e_E(1))/k/T)'/Zvibr(1);
            M.s_e(2)*exp(-(M.ev_i{2}+M.ev_0(2)+M.e_E(2))/k/T)'/Zvibr(2);
            M.s_e(3)*exp(-(M.ev_i{3}+M.ev_0(3)+M.e_E(3))/k/T)'/Zvibr(3)];
out=n;
end