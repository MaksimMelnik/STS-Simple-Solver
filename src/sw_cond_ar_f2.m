% программа пересчета параметров на УВ 
% смесь Ar-CO

% 18.11.2019
% 18.07.2022 исправленно
function out=sw_cond_ar_f2(T, v0, p0, f, num_vibr_levels)
k = 1.380649e-23;
syms p v T0 % dim-less variables
mCO=4.651236272601599e-26;
mAr=6.633521356992e-26;
F=mAr*(1-f)+mCO*f; 

S = [p*v/T == p0*v0/T0,...
     F*p/(k*T)*v^2+p==F*p0/(k*T0)*v0^2+p0,...
     2.5*k*T + f*k*T + Evibr(T0,num_vibr_levels)*f + F*v^2*0.5==...
        2.5*k*T0 + f*k*T0 + Evibr(T0,num_vibr_levels)*f + F*v0^2*0.5];

N = vpasolve(S,[p,v,T0],[200000,700,300]);
X = [N.p, N.v, N.T0];
out = double(X);
end