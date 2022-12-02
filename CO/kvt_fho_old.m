function kvt = kvt_fho_old(t, M1, M2, ind_e)
% программа расчета коэффициентов скорости VT обмена
% соотношения осредненных к.с. Adamovich et al. 1998
% внесены поправки согласно коду fho.f
% AB(i) + C = AB(f) + C
% 13/06/2017
% Модификация для универсальности и без глобальных переменных. Только для
% одноквантовых переходов (заменены факториалы на константы).
% t is temperature (K), M1 is the molecule under consideration, e_i is
% vibrational energy (J), num_v_l is the number of vibrational states, 
% ind_p is the indicator of a parner-particle, ind_e is the electronic
% state of M1.
% 02.12.2022

% константы
h = 6.6261*10^(-34); % постоянная Планка (Дж*сек)
h_bar = h/(2*pi);
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
N_a = 6.0221e23; % постоянная Авогадро (1/моль)
c = 299792458; 

alpha = 4e10;
depth = 200.0;
e_i=M1.ev_i{ind_e};
% газокинетический радиус упругого столкновения, CO-CO
% r0 = 3.65e-10; % m, для CO
r0=(M1.diameter+M2.diameter);
% сечение столкновения
sigma = pi*r0^2; % модель твердых сфер
% reduced collision mass
mu = M1.m_mass*M2.m_mass./(M1.m_mass+M2.m_mass) * 1e-3 / N_a;
% частота столкновений
Z = sigma*sqrt(8*k*t./(pi*mu)) * 1e6; % -> cm^3/sec

% mpar = ma*mc/mb./(ma+mb+mc);
% Svt = 2*mpar./(1+mpar).^2; % Svt = 1/pi;
Svt=0.5;
% nu = 1; % столкновение молекула-молекула
% nu = 0; % столкновение молекула-атом
nu=1;
if M2.fr_deg_c==3
    nu=0;
end

% l=68;
% Lmax=num_v_l-1;
Lmax=M1.num_vibr_levels(ind_e)-1;
i=2:M1.num_vibr_levels(ind_e);

f = i-1;

delE = (e_i(i)-e_i(f))/h/c*1.4388e-2;
theta = abs(delE./(i-f));               % K
om = theta*k/h_bar;                     % sec^-1
theta1 = 4*pi^2*om.^2.*mu/alpha^2/k;    % K

s = abs(i-f);
sf=s;
ns=i;

% эффективная скорость ЛТ
vm0 = (2.*pi*om.*s*k*t/alpha./mu).^(1/3);
const = 1./s .* (nu+2*ns.^(1./s)./(s+1)).*Svt.*theta1./theta;
const2 = (ns.^(1./s)./(s+1).*Svt.*theta1./theta).^2;
x0 = 2;
Cvt=zeros(1,Lmax);
% for jj=1:Lmax
%    fun = @(x) x-nthroot(1-const(jj).*exp(-2*pi.*om(jj)./...
%        (alpha*vm0(1,jj)*x))-const2(jj).*exp(-4*pi*om(jj)./...
%        (alpha*vm0(1,jj)*x)), 3);
%    Cvt(1,jj)= fzero(fun,x0);
% end
for jj=1:Lmax
   fun = @(x) x.^3-(1-const(jj).*exp(-2*pi.*om(jj)./...
       (alpha*vm0(1,jj)*x))-const2(jj).*exp(-4*pi*om(jj)./...
       (alpha*vm0(1,jj)*x)));
   Cvt(1,jj)= fzero(fun,x0);
end

delta = (1-Cvt.^3)./Cvt.^3 * 2*pi.*om/alpha./vm0./Cvt;
    phi = 2/pi * atan(sqrt(2*k*depth/mu)./vm0(1,:));


rate = ns.*sqrt(2*pi./(3+delta)).*s.^(1/3)./sf.^2;
    rate(1,:) = rate(1,:) .* Cvt(1,:) .* (Svt.*theta1(1,:)./theta).^s...
                                                 .*(theta1(1,:)/t).^(1/6);
rate = rate.*exp(- s.^(2/3).*(theta1/t).^(1/3).*(0.5.*Cvt.^2+1./Cvt).*...
                                            (1-phi).^(2/3)-s.*(1-Cvt.^3));
rate = rate .* exp(theta.*s/2/t);
    kvt(1,:) = rate(1,:) * Z;
kvt=kvt/1e6;

