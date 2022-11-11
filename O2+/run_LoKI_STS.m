function run_LoKI_STS(T, n1, nO2a, nO2b, ne, nO, nO1D, n0, EdN, t)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
nAll=sum(n1)+nO2a+nO2b+nO+nO1D;
oxy_saver_STS(nAll*n0*k*T, T, ne*n0, n1*n0, ...
    nO2a*n0, nO2b*n0, nO*n0, nO1D*n0, EdN);
cd LoKI-master\Code
data4me=lokibcl_f('Oxygen\oxygen.in');
cd ..\..
disp([num2str(t) ' run LoKI-B'])
end