function out = R_diss_old(T, ni_b, nm_b, nac_b, nao_b, naAr_b, modD,...
    ind_Arr, n0)
num_vibr_levels=69;    
N_a=6.02214076e23; V_C = 299792458; V_K = 1.380649e-23;
V_H = 6.626070041e-34; V_CO_diss_E =  1.77716869535000e-18;
V_CO_mass = 4.651236272601599e-26;
V_C_mass=1.994473440944e-26; V_O_mass = 2.6567628316576e-26;
V_CO_form_E=-1.86384779872000e-19; V_C_form_E = 1.18095841117000e-18;
    V_O_form_E = 4.098045681049634e-19;
if nargin<7
   ind_Arr=1;
end
%ArrA = 3.819239702291926e-10; ArrN = -1;   % Park for all
switch ind_Arr
    case 1
        ArrACO=2.3e20/N_a*1e-6; ArrAC=3.4e20/N_a*1e-6; 
        ArrAO=3.4e20/N_a*1e-6; 
        ArrNCO=-1; ArrNC=-1; ArrNO=-1; 
        ArrAAr=2.3e19/N_a*1e-6; ArrNAr=-1;      % Park
%         ArrAAr=4.3e27/N_a*1e-6; ArrNAr=-3.1;    % Mick
    case 2
        ArrACO=2.8e21/N_a*1e-6; ArrAC=1.4e22/N_a*1e-6; 
        ArrAO=2.1e22/N_a*1e-6; 
        ArrNCO=-1.39; ArrNC=-1.39; ArrNO=-1.39;
        ArrAAr=1.4e21/N_a*1e-6; ArrNAr=-1.39;
end

kd_eq = exp(-V_CO_diss_E /(V_K*T));% m^3/sec
e_i=(levels_e(num_vibr_levels))';
e_i=e_i-e_i(1);
ZvT = sum(exp(-e_i/(V_K*T)));
% parameter of TM model
switch modD
    case 2
        U = V_CO_diss_E/(V_K*6);
    case 3
        U = 3*T;
    case 4
        U = Inf;
    case 5
        U = V_CO_diss_E/V_K *(0.15+T/20000);
end
ZvU = sum(exp(e_i/(V_K*U)));
% non-equilibrium factor
Z = ZvT / ZvU * exp(e_i/V_K*(1/T + 1/U));
% dis. rates
kd = kd_eq' * Z'; % m^3/sec
sigma = 1;
Be=193.128087;
Theta_r = Be*V_H*V_C/V_K;
Z_rot = T./(sigma.*Theta_r);
Kdr2=(V_CO_mass/V_C_mass/V_O_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2)...
    *Z_rot *exp(-(e_i'+(V_CO_form_E-V_C_form_E-V_O_form_E))/V_K/T);

kr= kd .* Kdr2 * n0;
kd_CO=kd .*ArrACO.*T.^ArrNCO;
kd_C=0;%kd.*ArrAC.*T.^ArrNC;
kd_O=0;%kd.*ArrAO.*T.^ArrNO;
kd_Ar=0;%kd.*ArrAAr.*T.^ArrNAr;
kr_CO=kr.*ArrACO.*T.^ArrNCO;
kr_C=0;%kr.*ArrAC.*T.^ArrNC;
kr_O=0;%kr.*ArrAO.*T.^ArrNO;
kr_Ar=0;%kr.*ArrAAr.*T.^ArrNAr;

RD = nm_b * (nac_b*nao_b*kr_CO-(ni_b').*kd_CO) + ...
    nac_b * (nac_b*nao_b*kr_C-(ni_b').*kd_C)+ ...
    nao_b * (nac_b*nao_b*kr_O-(ni_b').*kd_O)+ ...
    naAr_b* (nac_b*nao_b*kr_Ar-(ni_b').*kd_Ar);
out = RD';
gg=-(e_i'+(V_CO_form_E-V_C_form_E-V_O_form_E));
disp(gg(1))
% pause
end