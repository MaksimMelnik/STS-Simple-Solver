function kvv = kvv_fho(t, O2)
%i1->i1-1
%i2->i2+1
% global O2
l=O2.num_vibr_levels(1);
Lmax=l-1;
% константы
h = 6.6261*10^(-34); % постоянная Планка (Дж*сек)
hbar = h/(2*pi);
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
N_a = 6.0221e23; % постоянная Авогадро (1/моль)
c = 2.99*10^10; % скорость света (м/сек)

T = t;


% % N2-N2
% ome = 3394.2; % K
% omexe = 20.612; %K
% ma = 14;
% mb = 14;
% mc = 28; % N2
% alpha = 4e10; % m^-1
% depth = 200; % K

% спектроскопические постоянные O2, K
ome = 2273.7;
omexe = 17.369;
% молекулярный вес атомов в O2 и партнера, г/моль
ma = 16.0;
mb = 16.0;
mc = 32.0;
alpha = 4e10;
depth = 200.0;

gaskin=3.458e-8^2;
mu = (ma+mb)*mc/(ma+mb+mc) * 1e-3 / N_a; % kg
mpar = ma*mc/mb/(ma+mb+mc);
d8pi = sqrt(8*pi);
velt = sqrt(k*T/mu); 
delv = 0.01*velt;
freq = gaskin*100*d8pi*velt;
fvpar = d8pi/velt;
 
i1=1:Lmax;
f1=i1-1;
i2=f1;
f2=i1;

% Evibi1 = (ome*(i1+0.5) - omexe*(i1+0.5).^2)';
% Evibi2 = ome*(i2+0.5) - omexe*(i2+0.5).^2;
% Evibf1 = (ome*(f1+0.5) - omexe*(f1+0.5).^2)';
% Evibf2 = ome*(f2+0.5) - omexe*(f2+0.5).^2;

Evibi1=(O2.e_i(1, 2:end)'+O2.e_0)*1.4388e-2/h/c*1e2;
Evibi2=(O2.e_i(1, 1:end-1)+O2.e_0)*1.4388e-2/h/c*1e2;
Evibf1=Evibi2';
Evibf2=Evibi1';

s = abs(i1-f1);
sf = factorial(s);
% ns1 = (factorial(max(i1,f1))./factorial(min(i1,f1)));
% ns2 = (factorial(max(i2,f2))./factorial(min(i2,f2)));
% disp('ns')
ns1=1:41;
ns2=1:41;
    eee1 = abs((Evibi1-Evibf1)./(i1-f1)');
    eee2 = abs((Evibi2-Evibf2)./(i2-f2));
    Evib = (eee1+eee2)/2;
delE = Evibi1-Evibf1+Evibi2-Evibf2;
delom = 2*pi*delE*c/1.439;
omeg = 2*pi*Evib*c/1.439;

% ro2 = alpha^2*k*T/(2*omeg^2*mu)/27;
ro2 = alpha^2*k*T./(2*omeg.^2*mu)/16;       
 
% Resonance defect correction (Keck, Carrier, Adamovich)
popc = 4*pi^2 * mu * omeg.^2 / alpha^2 / k;
popc = sqrt(popc/T) .* abs(delE)./Evib / sqrt(8);
popc = (2/3)*popc./s;
popc = exp(-2/3 * popc);
popc = 0.5*(3-popc).*popc;

% FHO rates for the exact resonance with the correction
zompl0 = (ns1.^(1./s))' .* ns2.^(1./s);
zompl1 = zompl0.*ro2;
zompl2 = zompl1.^s./sf;
zompl3 = (1+2*zompl1./(s+1)).^(s+1);
zompl4 = zompl2./zompl3;
kvv = zompl4 .* popc * freq;
%вот ето не факт, но у Адамовича так
kvv=kvv.*exp(0.5*delE/T)/1e6;
end