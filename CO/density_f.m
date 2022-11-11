function out=density_f(T, n1, num_vibr_levels)
V_K = 1.380649e-23;
e_i = levels_e(num_vibr_levels);
Zv=sum(exp(-(e_i-e_i(1))/T/V_K));
n=n1/Zv.*exp(-(e_i-e_i(1))/T/V_K);
out=n;
end