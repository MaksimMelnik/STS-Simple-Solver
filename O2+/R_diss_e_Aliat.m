function out = R_diss_e_Aliat(M, n_m, n_a1, n_a2, n_p, kd_eq, T, n0, ...
    ind_U)
% Универсальная функция расчёта диссоциации по Алиату для
% столкновения электронно-возбуждённых молекул с электронами, где 
% коэффициент скорости диссоциации считается вне.
% 24.11.2020
V_C = 299792458; V_K = 1.380649e-23; V_H = 6.626070041e-34;
if nargin<9
   ind_U=2;
end

% parameter of TM model
switch ind_U
    case 2
        U = M.diss_e/(V_K*6);
    case 3
        U = 3*T;
    case 4
        U = Inf;
    case 5
        U = M.diss_e/V_K *(0.15+T/20000);
end
    % Aliat
Zel=M.s_e*exp(-M.e_E'/(V_K*T));
divsum=0;
V=[];
Kdr=[];
Theta_r = M.Be*V_H*V_C/V_K;
Z_rot = T./(M.sigma.*Theta_r);
for i=1:M.num_elex_levels
    ZvibT=sum(exp(-M.e_i(i, 1:M.num_vibr_levels(i))/(V_K*T)));
    ZvibU=sum(exp(M.e_i(i, 1:M.num_vibr_levels(i))/(V_K*U)));
    divsum=divsum+M.s_e(i)*exp(M.e_E(i)/(V_K*U))*ZvibU/ZvibT;
    V=[V exp((M.e_i(i, 1:M.num_vibr_levels(i))+M.e_E(i))/V_K*(1/T+1/U))];
    Kdr=[Kdr Z_rot(i)*M.s_e(i)*...
        exp(-(M.e_i(i,1:M.num_vibr_levels(i))+M.e_E(i)-M.diss_e)/V_K/T)];
end
V=Zel*V/divsum;
Kdr=Kdr*(M.mass/M.mltpl_atoms_mass)^(3/2)*V_H^3*(2*pi*V_K*T)^(-3/2);
% kd_new=V*kd_eq;

kd=kd_eq;
% kd=kd_new;
Kdr2=0; %Kdr;
kr= kd .* Kdr2 * n0;
kr2=kd .* Kdr *  n0;
RD = n_p * (n_a1*n_a2*kr-(n_m').*kd);
% disp('inside')        % try to fix recombination
% disp(['n_a ' num2str(n_a1) ' na*kr2 ' num2str(n_a2*kr2) ' kr ' num2str(kr2) ' kd ' num2str(kd)])
% disp('Kdr')
% disp(Z_rot)
% i=1;
% disp(exp(-(M.e_i(i,1:M.num_vibr_levels(i))+M.e_E(i)-M.diss_e)/V_K/T))
out = RD';
end