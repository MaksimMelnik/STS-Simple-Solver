function tau=draw_tau_VT_Ar(temp, v0)
% рисуем график колебательной релаксации по полученным данным
% 01.04.2020
V_K = 1.380649e-23; 
num_vibr_levels = 69;%41;
    gg=(temp(:,end-2)-temp(:,end-1))./temp(:,end-2);
    ii=find(abs(gg)<1e-3);
    tau=temp(ii(1),1)/v0*1e6;
    ptau=tau*sum(temp(1,2:num_vibr_levels+4),2)*V_K*temp(1,end-2)/101325;
    Ttau=temp(1,end-2);
    figure(3)
    semilogy(Ttau, ptau/1e6, 'p')%, 'MarkerSize', '12', ...
                                    % 'MarkerFaceColor', 'k')
end