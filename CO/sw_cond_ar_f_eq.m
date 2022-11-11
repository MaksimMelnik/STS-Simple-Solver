% программа пересчета параметров на УВ для смести Ar/CO/C/O по давлению 
% перед УВ, скорости УВ, начальной концентрации CO и температуре в конце
% релаксации
% 18.04.2020
function out=sw_cond_ar_f_eq(T, v0, p0, f, num_vibr_levels)
k = 1.380649e-23;
h = 6.626070041e-34;
N_a=6.02214076e23;
syms nco nc na v T0 % dim-less variables
mCO=4.651236272601599e-26;
mAr=6.633521356992e-26;
mC=1.994473440944e-26;
mO=2.6567628316576e-26;
V_CO_diss_E = 1.77716869535000e-18; V_CO_form_E=-1.86384779872000e-19;
V_C_form_E = 1.18095841117000e-18; V_O_form_E = 4.098045681049634e-19;
ICO=1.456038183037318e-46; % момент инерции
F=mAr*(1-f)+mCO*f;

ArrACO=2.3e20/N_a*1e-6; ArrAC=3.4e20/N_a*1e-6; 
        ArrAO=3.4e20/N_a*1e-6; 
        ArrNCO=-1; ArrNC=-1; ArrNO=-1; 
        ArrAAr=2.3e19/N_a*1e-6; ArrNAr=-1;
kdCO=ArrACO*T^ArrNCO*exp(-V_CO_diss_E /(k*T));
kdC=ArrAC*T^ArrNC*exp(-V_CO_diss_E /(k*T));
kdO=ArrAO*T^ArrNO*exp(-V_CO_diss_E /(k*T));
kdAr=ArrAAr*T^ArrNAr*exp(-V_CO_diss_E /(k*T));
e_i=levels_e(num_vibr_levels);
Zvibr=sum(exp(-e_i/k/T));
Zrot=8*pi^2*ICO*k*T/h^2;
Zint=Zrot*Zvibr;
db=(mCO/mC/mO)^1.5*h^3/(2*pi*k*T)^1.5*Zint*exp(V_CO_diss_E/(k*T));
krCO=db*kdCO;
krC=db*kdC;
krO=db*kdO;
krAr=db*kdAr;

na=(nco+nc)/f*(1-f);

S = [(nco+2*nc+na)*k*v == p0*v0/T0,...
     (F/(k*T)*v^2+1)*(nco+2*nc+na)*k*T==F*p0/(k*T0)*v0^2+p0,...
     ...2.5*k*T + f*k*T + Evibr(T,num_vibr_levels)*f + F*v^2*0.5==...
        ...2.5*k*T0 + f*k*T0 + Evibr(T0,num_vibr_levels)*f +F*v0^2*0.5,...
     (2.5*(nco+2*nc+na)*k*T+nco*k*T+Evibr(T,num_vibr_levels)*nco+...
        nco*V_CO_form_E+nc*(V_C_form_E+V_O_form_E))/...
        (mCO*nco+nc*(mC+mO)+na*mAr)+v^2*0.5==...
        (2.5*k*T0 + f*k*T0 + Evibr(T0,num_vibr_levels)*f + ...
        f*V_CO_form_E)/F+v0^2*0.5,...
     nco*(nc^2*krCO-nco*kdCO)+na*(nc^2*krAr-nco*kdAr)+... % тут R==0
        nc*(nc^2*(krC+krO)-nco*(kdC+kdO))==0 ...          % и тут
     ];

N = vpasolve(S,[nco,nc,v,T0],[1e21,1e21,700,600]);
X = [N.nco, N.nc, N.v, N.T0];
out = double(X);
end