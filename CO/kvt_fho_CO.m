function kvt = kvt_fho_CO(t, e_i, num_v_l, ind_p)
% программа расчета коэффициентов скорости VT обмена
% соотношения осредненных к.с. Adamovich et al. 1998
% внесены поправки согласно коду fho.f
% AB(i) + C = AB(f) + C
% 13/06/2017
% модификация для CO и кода без глобальных переменных. Только для
% одноквантовых переходов (заменены факториалы на константы).
% Создана неизвестно когда. доработка продолжилась
% 15.08.2021

% константы
h = 6.6261*10^(-34); % постоянная Планка (Дж*сек)
h_bar = h/(2*pi);
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
N_a = 6.0221e23; % постоянная Авогадро (1/моль)
c = 299792458; 

alpha = 4e10;
depth = 200.0;
 ma=12;         % C
 mb=16;         % O
 mAr=39.95;     % Ar
 mc(1)=28;      % CO
 mc(2)=ma;
 mc(3)=mb;
 mc(4)=mAr;
 mc(5)=ma*2;    % C2
% газокинетический радиус упругого столкновения, CO-CO
r0 = 3.65e-10; % m, для CO
% сечение столкновения
sigma = pi*r0^2; % модель твердых сфер
% приведенная масса CO-CO/C/O
mu = (ma+mb)*mc(ind_p)./(ma+mb+mc(ind_p)) * 1e-3 / N_a;
% частота столкновений
Z = sigma*sqrt(8*k*t./(pi*mu)) * 1e6; % -> cm^3/sec

% mpar = ma*mc/mb./(ma+mb+mc);
% Svt = 2*mpar./(1+mpar).^2; % Svt = 1/pi;
Svt=0.44;
% nu = 1; % столкновение молекула-молекула
% nu = 0; % столкновение молекула-атом
nu=[1 0 0 0 1];

% l=68;
Lmax=num_v_l-1;
i=2:num_v_l;

f = i-1;

delE = (e_i(i)-e_i(f))/h/c*1.4388e-2;
theta = abs(delE./(i-f)); % K
om = theta*k/h_bar; % sec^-1
% theta1 = 4*pi^2*om.^2.*mu/alpha^2/k; % K
theta1(1,:) = 4*pi^2*om.^2.*mu(1)/alpha^2/k; % K
% theta1(2,:) = 4*pi^2*om.^2.*mu(2)/alpha^2/k; % K
% theta1(3,:) = 4*pi^2*om.^2.*mu(3)/alpha^2/k; % K

s = abs(i-f);
% sf = factorial(s);
sf=s;
% ns = (factorial(max(i,f))./factorial(min(i,f)));
ns=i;

% эффективная скорость ЛТ
vm0(1,:) = (2.*pi*om.*s*k*t/alpha./mu(1)).^(1/3);
% vm0(2,:) = (2.*pi*om.*s*k*t/alpha./mu(2)).^(1/3);
% vm0(3,:) = (2.*pi*om.*s*k*t/alpha./mu(3)).^(1/3);
%
% for j=1:3
% const(j,:) = 1./s .* (nu(j)+2*ns.^(1./s)./(s+1)).*Svt.*theta1(j,:)./theta;
% const2(j,:) = (ns.^(1./s)./(s+1).*Svt.*theta1(j,:)./theta).^2;
% end
const = 1./s .* (nu(ind_p)+2*ns.^(1./s)./(s+1)).*Svt.*theta1./theta;
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

% gg=max(abs(Cvt-Cvt3)./Cvt);
% if gg>0
% disp(gg)
% end

%    fun2 = @(x) x.^3-(1-const.*exp(-2*pi.*om./...
%        (alpha*vm0.*x))-const2.*exp(-4*pi*om./...
%        (alpha*vm0.*x)));
% options=optimoptions('fsolve','Display','none');
%    Cvt= fsolve(fun2, const*0+x0, options);
   
%    x=Cvt;
%    temp=1-const.*exp(-2*pi.*om./...
%        (alpha*vm0.*x))-const2.*exp(-4*pi*om./...
%        (alpha*vm0.*x));
%    disp([min(temp) max(temp)])

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

