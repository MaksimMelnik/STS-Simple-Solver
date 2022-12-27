%% calculate incubation times
arr_i_ini=[2 4 5];
for i_ini=arr_i_ini
 for i_dis=1:3
  for i_U=2
%       disp([i_ini i_dis])
x=res_all(i_ini, i_dis, i_U).temp(:, end); % time, musec
            y=res_all(i_ini, i_dis, i_U).temp(:,1+num_C); % nC
            d1=diff(y)./diff(x);
            x1=(x(1:end-1)+x(2:end))/2;
            ii=find(d1>max(d1)*0.85);
            ii=ii(x1(ii)<50);
            i1=ii(1);
            i2=ii(end);
            x0=(x1(i2)+x1(i1))/2;
            interY=spline(x(i1:i2), y(i1:i2), x0);
            interd1=spline(x1(i1:i2), d1(i1:i2), x0);
            tau=(-interY/interd1+x0); % собственно tau, в мкс, если 
                                         % правильно понял
tau_arr(i_ini, i_dis, i_U)=tau;
% figure  % рисуем касательную и точку инкубации
%     plot(x, y)
%     hold on
%     dx0=30;
%     kas=[interd1*(-dx0)+interY interY interY+interd1*dx0];
%     plot([x0-dx0 x0 x0+dx0],kas)
%     plot([-10 10], [0 0], 'k')
%     plot(tau, 0, '*')
%     xlim([0 20])
% pause
  end
 end
 T(i_ini)=res_all(i_ini, 1, 2).temp(1,1+num_T);
end

semilogy(T, tau_arr(:, 1, 2))
xlim([7000 16000])
ylim([1e-2 1e2])