function k_down = kvv_billing(t, M1, M2, ind_e1, ind_e2, anharm)
% kdown - array of k_(i->i-1)^(j->j+1)
l=M1.num_vibr_levels(ind_e1);           % number of vibr. lvls
Lmax = l - 1;
if M1.name() == "N2" && M2.name() == "N2"
    k_down = zeros(Lmax, Lmax);
    d_vv = 6.8/sqrt(t);
    j = 0:Lmax-1;
    for i = 1:Lmax
        k_down(i, :) = (j+1).*2.5*(10^(-14))*i*((t/300)^(3/2))...
            .*exp((i-1-j)*d_vv).*((3/2)-(1/2)*exp(-(abs(i-1-j))*d_vv)); % i-1<=k?
    end
elseif M1.name() == "O2" && M2.name() == "O2"
    k_down = zeros(l-1, l-1);
    d_vv = 2.4/sqrt(t);
    j = 0:Lmax-1;
    for i = 1:Lmax
        k_down(i, :) = 2.8*(10^(-18))*i.*(j+1)*((t)^(3/2)).*exp((i-1-j)*d_vv);
    end
else
    k_down = kvv_ssh(t, M1, M2, ind_e1, ind_e2, anharm);
end
k_down = k_down * 1e-6;
end
