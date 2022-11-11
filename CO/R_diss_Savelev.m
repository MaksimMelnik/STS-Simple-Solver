function out = R_diss_Savelev(M, n_m, n_a1, n_a2, n_p, Coll, T, ...
                                            n0, ind_Arr, ind_U, ind_exc)
% Универсальная функция расчёта диссоциации по Савельеву для
% электронно-возбуждённых молекул (статья с Поляховских).
% 08.09.2021
c = 299792458; k = 1.380649e-23; h = 6.626070041e-34;
if nargin<9
   ind_Arr=1;
   ind_U=2;
   ind_exc=1;
end

% В переменной Coll хранится информация о столкновении. ArrA и ArrN --
% параметры в законе Аррениуса. ArrA(1) -- по Парку, 2 -- по Ибрагимовой,
% 3 -- по МакКензи.
kd_eq = Coll.ArrA(ind_Arr) * T^Coll.ArrN(ind_Arr)*...
                            exp(-(M.diss_e + M.e_E)/(k*T)); % m^3/sec
if ind_Arr==4
    kd_eq = kd_eq*0 + Coll.ArrA(ind_Arr);
end

    % parameter of TM model
U=0;%M.diss_e*0;
switch ind_U
    case 2
        U = (M.diss_e+M.e_E)/(k*6);
    case 3
        U = U+3*T;
    case 4
        U = U+Inf;
%     case 5
%         U = U+M.diss_e/k *(0.15+T/20000);
end
% dE = activation + product - reagent
% dE= Ea + Ep - Er
    % упаковываем энергию всех электронно-колебательных уровней в один
    %   массив
Eli=[   M.ev_i{1}'+M.e_E(1)+M.ev_0(1);           % эл.-кол. энергия для X
        M.ev_i{2}'+M.e_E(2)+M.ev_0(2);...          для a3П
        M.ev_i{3}'+M.e_E(3)+M.ev_0(3)];          % для A1П
    % El - это энергия одного из реагентов - энергия активации
E1=[M.ev_i{1}'+M.e_E(1)+M.ev_0(1)-M.e_E(1)-M.diss_e(1);
    M.ev_i{2}'+M.e_E(2)+M.ev_0(2)-M.e_E(2)-M.diss_e(2);
    M.ev_i{3}'+M.e_E(3)+M.ev_0(3)-M.e_E(3)-M.diss_e(3)];
    % энергия 1 реагента + энергия 2 реагента - энергия активации,
    %   т.е. Er-Ea.
    % Этот кусок немного мудрённый, но так проще записывать.
Er_a=Eli'+E1;
    % Обезразмеренное равновесное распределение по Больцману
    %   по электронно-колебательным состояниям.
neq=density_f_exc(T, 1, M);
    % Не важно, сложности с размерностями для U=D/6k.
if ind_U==2
 U=[    M.ev_i{1}'*0 + M.e_E(1)+M.diss_e(1); ...
        M.ev_i{2}'*0 + M.e_E(2)+M.diss_e(2);...
        M.ev_i{3}'*0 + M.e_E(3)+M.diss_e(3)]/(k*6);
end
    % Не важно, сложности с размерностью.
kd_eq_B=[   M.ev_i{1}'*0+kd_eq(1); ...
            M.ev_i{2}'*0+kd_eq(2);...
            M.ev_i{3}'*0+kd_eq(3)];
s_e=[   M.ev_i{1}'*0+M.s_e(1); ...
        M.ev_i{2}'*0+M.s_e(2);...
        M.ev_i{3}'*0+M.s_e(3)];
Theta_r = M.Be*h*c/k;
Z_rot = T./(M.sigma.*Theta_r);
Z_rot_v=[   M.ev_i{1}'*0+Z_rot(1); ...
            M.ev_i{2}'*0+Z_rot(2);...
            M.ev_i{3}'*0+Z_rot(3)];
switch ind_exc
 case 1
  dE=-E1;
  dEexp=exp(-dE.*heaviside(dE)/k.*(1/T+1./U));
  B=kd_eq_B/sum(dEexp.*neq);
  kd=B.*dEexp;
  Kdr=s_e/M.mltpl_atoms_s_e*(M.mass/M.mltpl_atoms_mass)^1.5...
    *h^3*(2*pi*k*T)^(-1.5).*Z_rot_v.*exp(-(E1)/k/T);
  kr=0;%kd .* Kdr * n0;
  RD = kr*n_a1*n_a2*sum(n_m) -kd.*n_m*sum(n_m);
 case 2
  dE=M.ev_0(1)-Er_a;
  dEexp=exp(-dE.*heaviside(dE)/k.*(1/T+1./U));
  B=kd_eq_B/sum(dEexp.*neq.*neq', 'all');
  kd=B.*dEexp;
  Kdr=s_e/M.mltpl_atoms_s_e...%.*s_e'./s_e(1)...
    *(M.mass/M.mltpl_atoms_mass)^1.5*h^3*(2*pi*k*T)^(-1.5)...
    .*Z_rot_v...%.*Z_rot_v'./Z_rot_v(1)...
    .*exp(-(E1)/k/T);%.*exp((Eli'-M.e_0(1))/k/T);
  kr=kd .* Kdr * n0;
%   kr=0;
  RD = kr*n_a1*n_a2*n_m(1) -kd.*n_m.*n_m';
 case 3
    % Энергия продукта. ' и reshape подгоняют под нужные размерности.
  Ep=reshape(Eli, 1, 1, []);
    % Ep-Ee+Ea
  dE=Ep-Er_a;
  dEexp=exp(-dE.*heaviside(dE)/k.*(1/T+1./U));
  B=kd_eq_B/sum(dEexp.*neq.*neq', 'all');
  kd=B.*dEexp;
  Kdr=s_e/M.mltpl_atoms_s_e.*s_e'./reshape(s_e, 1, 1, [])...
    *(M.mass/M.mltpl_atoms_mass)^1.5*h^3*(2*pi*k*T)^(-1.5)...
    .*Z_rot_v.*Z_rot_v'./reshape(Z_rot_v, 1, 1, [])...
    .*exp(-(E1)/k/T).*exp(-(Eli'-reshape(Eli, 1, 1, []))/k/T);
  kr=kd .* Kdr * n0;
%   kr=0;
  RD = kr*n_a1*n_a2.*reshape(n_m, 1, 1, []) -kd.*n_m.*n_m';
end
out=RD;
end