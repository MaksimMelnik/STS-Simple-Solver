function k_down = kvt_billing(t, M1, M2, ind_e, anharm)

l=M1.num_vibr_levels(ind_e);    % number of vibr. lvls
Lmax = l-1; % number of maximum level
i = (1:Lmax);   
if M1.name() == "N2"
 if M2.name() == "N2"
  k10 = exp(-3.24093 - 140.69597/t^0.2);
  d_vt = 0.26679 - 6.99237 * 10 ^ (-5) * t + 4.70073 * 10 ^ (-9) * t ^ 2;
  k_down = i .* k10 .* exp((i-1) * d_vt);
 elseif  M2.name() == "N"
  b0 = -25.708 - 5633.1543/t;
  b1 = -0.1554 + 111.3426/t;
  b2 = 0.0054 - 2.189/t;
  c0 = 0.0536 + 122.4835/t;
  c1 = 0.0013 - 4.2365/t;
  c2 = -1.197 * 10^(-4) + 0.0807/t;
  k_down = zeros(size(i));
  for j = 1:Lmax
   %i1 = j-1;
   %k_down(j) = exp(b0+b1*(j-i1)+b2*(j-i1)^2+j*(c0+c1*(j-i1)+c2*(j-i1)^2));
   k_down(j) = exp(b0+b1+b2+j*(c0+c1+c2));
  end
 elseif M2.name() == "O2"
     %N2.num_elex_levels=1; ?
     k_down = kvt_billing(t, M1, N2, ind_e, anharm);
 elseif M2.name() == "O"
     k_down = kvt_billing(t, M1, N, ind_e, anharm);
 elseif M2.name() == "NO"
     k_down = zeros(Lmax);
 end
elseif M1.name() == "O2"
 if M2.name() == "O2"
  k10 = t/(1.8*10^12*exp(122/t^(1/3))*(1-exp(-2273.7/t)));
  d_vt = 2.99/sqrt(t);
  k_down = i.*k10.*exp((i-1) * d_vt);
 elseif M2.name() == "O" % слишком маленькие коэффициенты
  k_down = zeros(size(i));
  for j = 1:Lmax
   k_down(j) = (3*j-2)*(M1.ev_i{ind_e}(2))*(7*10^(-14))/(0.124*exp(30/t^(1/3))); %M1.ev_i{ind_e}(2)?
  end
 elseif M2.name() == "N2"
     k_down = kvt_billing(t, M1, O2, ind_e, anharm);
 elseif M2.name() == "N"
     k_down = kvt_billing(t, M1, O, ind_e, anharm);
 elseif M2.name() == "NO"
     k_down = zeros(Lmax);
 end
end
k_down = k_down/1e6;
end

