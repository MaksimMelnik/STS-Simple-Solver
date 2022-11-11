function kvt = kvt_fho(t, O2)
% программа расчета коэффициентов скорости VT обмена
% соотношения осредненных к.с. Adamovich et al. 1998
% внесены поправки согласно коду fho.f
% AB(i) + C = AB(f) + C
% 13/06/2017

% константы
c = 2.99*10^8; % скорость света (м/сек)
h = 6.6261*10^(-34); % постоянная Планка (Дж*сек)
h_bar = h/(2*pi);
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
% N_a = 6.0221e23; % постоянная Авогадро (1/моль)
mu=[2.656880490194451e-26 1.771253660129634e-26];
mpar=[5.000000000000000e-01 3.333333333333333e-01];
alpha =     4.000000000000000e+10;
depth =   200;
l=O2.num_vibr_levels(1);

% спектроскопические постоянные O2, K
ome = 2273.7;
omexe = 17.369;
% молекулярный вес атомов в O2 и партнера, г/моль
ma = 16.0;
mb = 16.0;
mc(1) = 32.0;
mc(2)=16.0;
%
% alpha = 4e10;
% depth = 200.0;
%  ma=12;
%  mb=16;
%  mc(1)=28;
% газокинетический радиус упругого столкновения, O2-O2
r0 = 3.458e-10; % m
% сечение столкновения
sigma = pi*r0^2; % модель твердых сфер
% приведенная масса O2-O2
% mu = (ma+mb)*mc./(ma+mb+mc) * 1e-3 / N_a;
% частота столкновений
Z = sigma*sqrt(8*k*t./(pi*mu)) * 1e6; % -> cm^3/sec

% mpar = ma*mc/mb./(ma+mb+mc);
Svt = 2*mpar./(1+mpar).^2; % Svt = 1/pi;
% nu = 1; % столкновение молекула-молекула
% nu = 0; % столкновение молекула-атом
nu=[1 0];

Lmax=l-1;
i=1:Lmax;

f = i-1;
% Ei = ome*(i+0.5) - omexe*(i+0.5).^2;
% Ef = ome*(f+0.5) - omexe*(f+0.5).^2;
Ei=O2.e_i(1, 2:end)*1.4388e-2/h/c;         % e_i(0)=0
Ef=O2.e_i(1, 1:end-1)*1.4388e-2/h/c;

delE = Ei-Ef;
theta = abs(delE./(i-f)); % K
om = theta*k/h_bar; % sec^-1
% theta1 = 4*pi^2*om.^2.*mu/alpha^2/k; % K
theta1(1,:) = 4*pi^2*om.^2.*mu(1)/alpha^2/k; % K
theta1(2,:) = 4*pi^2*om.^2.*mu(2)/alpha^2/k; % K

s = abs(i-f);
sf = factorial(s);
ns = 1:Lmax;%(factorial(max(i,f))./factorial(min(i,f)));

% эффективная скорость ЛТ
vm0(1,:) = (2.*pi*om.*s*k*t/alpha./mu(1)).^(1/3);
vm0(2,:) = (2.*pi*om.*s*k*t/alpha./mu(2)).^(1/3);
%
for j=1:2
const(j,:) = 1./s .* (nu(j)+2*ns.^(1./s)./(s+1)).*Svt(j).*theta1(j,:)./theta;
const2(j,:) = (ns.^(1./s)./(s+1).*Svt(j).*theta1(j,:)./theta).^2;
end
% fun = @(x) x-nthroot(1-const*exp(-2*pi*om./(alpha*vm0*x)), 3);
x0 = 2;
Cvt=zeros(2,Lmax);
for j=1:2
   for jj=1:Lmax
       fun = @(x) x-nthroot(1-const(j,jj).*exp(-2*pi.*om(jj)./...
           (alpha*vm0(j,jj)*x))-const2(j,jj).*exp(-4*pi*om(jj)./...
           (alpha*vm0(j,jj)*x)), 3);
       Cvt(j,jj)= fzero(fun,x0);
   end
end

delta = (1-Cvt.^3)./Cvt.^3 * 2*pi.*om/alpha./vm0./Cvt;
for j=1:2
    phi(j,:) = 2/pi * atan(sqrt(2*k*depth/mu(j))./vm0(j,:));
end

rate = ns.*sqrt(2*pi./(3+delta)).*s.^(1/3)./sf.^2;
for j=1:2
    rate(j,:) = rate(j,:) .* Cvt(j,:) .* (Svt(j).*theta1(j,:)./theta).^s...
                                                  .*(theta1(j,:)/t).^(1/6);
end
rate = rate.*exp(- s.^(2/3).*(theta1/t).^(1/3).*(0.5.*Cvt.^2+1./Cvt).*...
                                             (1-phi).^(2/3)-s.*(1-Cvt.^3));
rate = rate .* exp(theta.*s/2/t);
for j=1:2
    kvt(j,:) = rate(j,:) * Z(j);
end
kvt=kvt/1e6;

