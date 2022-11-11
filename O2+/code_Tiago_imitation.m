% Код для расчёта 0d релаксации низкотемпературной плазмы в разряде.
% Пытаемся посчитать тест-кейс от Тиаго.
% На основе code_Tiago_test, но теперь стараемся велючить все процессы, 
% как у него в порядке влияния. Не важно, правильно ли мы подключаем)0
% потом разберёмся
% 22.11.2020

% ЗАДАЧИ:

% закон сохранения массы пофиксить (он с 1е-3 выполнялся)
    % сначала выровнять коэффициенты и сойтись по концентрациям с Тиаго
    % Попробовать пофиксить ошибку массы посчитав обратные по 
    %   детальному балансу
% разобраться с рекомбинацией: спросить у Тиаго про рек от е
% O(-,gnd) + O(3P) -> e + O2(X)     тут krate не совпадает в статье и
%   у Тиаго
% в некоторых процессах в статьях партнер -- О2(любой), а у Тиаго только
%   О2(Х). Поговорить с ЕВ, стоит ли вводить это.
% прочекать величины разные как в ЛоКИ чтобы были
% пофиксить вычисление Е для определения погрешности энергии
% попробовать считать ЛоКИ параллельно основному расчёту
% выяснить, как правильно высчитывать вклад диссоциации в изменение
%   температуры
% проверить пункт в ЛоКИ, который может оказаться потоком для тяжёлых
%   частиц, чтобы не суммировать Qin
% подключить Qin для разных процессов
% пересчитать e+O2(a)->e+2O(3P) по сечениям, а не как сейчас
% добавить Qin для процессов, не связанных с электронами
%   (а вроде португальцы не добавляют это)
    % добавить процессы
    % подключить
    % прогнать расчёт
    % спросить у ЕВ про диссоциацию (миро)
% подключить появление электронов от ионизации частиц
% добавить колебательную кинетику
% заново разобраться с диссоциацией по Алиату, почему там коэфф большие
%   вернее даже с рекомбинацией по нему
    % Физика А -- если оттуда статья Алиата с формулами, то у Алексея 
    %   спросить, какая там ошибка у Алиата
    % Вторая глава зеленой книжки, проверить, нужно ли учитывать 
    %   колебательные статсуммы в детальном балансе для поуровневого
    

fclose('all');
checkX=1e-20; %чтобы не слишком много выводить времени в комманд
ct=cputime-10;

%   Константы
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
R = 8.3144621;    % универсальная газовая постоянная
h = 6.626070041e-34;    % постоянная Планка (Дж/с)
HBAR = 1.054571800e-34;
c = 299792458;
Torr = 133.332;%.322368; % число Па в 1 Тор
alpha_FHO =4e10;
E_FHO =200;
SVT_FHO=0.44;

    % объявляем все параметры всем сортов газа и столкновений между ними
species_ini

            
          % p, Torr     fO2     fO2*    fO
incon_mas=[ 1           1       0       0
            1           0.7     0.1     0.2
            10          1       0       0
            10          0.7     0.1     0.2
        ];
ind=3;  % у третьего наибольший вклад от электронов
T=400;
p=incon_mas(ind, 1)*Torr;     % давление в паскалях. 1 или 10 Торр
EdN=50;                       % поле
n_all=p/(k*T);                % количество всех молекул
n1=n_all*incon_mas(ind, 2);   % количество молекул O2(X)
nO2a=n_all*incon_mas(ind, 3); % количество молекул O2(a)
nO=n_all*incon_mas(ind, 4);   % числовая плотность атомов кислорода O

ne=6.5e15;                  % числовая плотность электронов в м-3
Tv=T;                       % колебательная температура
is_e=1;
T0=T;                       % температура для обезразмеривания
n0=n1;                      % то же самое с n0 для обезразмеривания
tau=k/(p*pi*O2_O2_coll_diam*O2_O2_coll_diam)*sqrt(32*T/(6*R));

% % Распределение Больцмана
% 	Zv=sum(exp(-O2.e_i(1,:)/Tv/k));
%     n=n1/Zv.*exp(-O2.e_i(1,:)/Tv/k);
    n=n1;
    
E0=1.5*(n1+nO)*k*T+n1*k*T+(O2.e_i(1,:)+O2.e_0)*n'+O.form_e*nO;

num_O2=O2.num_vibr_levels(1);
y=n/n0;                         % сначала колебательные уровни
num_O2a=num_O2;
if O2.num_elex_levels>1
    num_O2a=num_O2+1;
    y(num_O2a)=n_all*incon_mas(ind, 3)/n0; % first electron excited state
end
num_O2b=num_O2a;
if O2.num_elex_levels>2
    num_O2b=num_O2a+1;
    y(num_O2b)=0/n0;                       % second electron excited state
end
num_O2p=num_O2b+1;
num_O=num_O2p+1;
num_O1D=num_O+1;
num_Op=num_O1D+1;
num_Om=num_Op+1;
num_O3=num_Om+1;
num_O3exc=num_O3+1;
num_e=num_O3exc+1;
num_T=num_e+1;
num_E_e=num_T+1;

y(num_O2p)=0;       % ion O2+
y(num_O)=nO/n0;     % атомы
y(num_O1D)=0;
y(num_Op)=0;        % ion O+
y(num_Om)=0;        % ion O-
y(num_O3)=0;
y(num_O3exc)=0;     % O3 excited
y(num_e)=ne/n0;     % электроны
y(num_T)=T/T0;
% y(num_E_e)=0;%Tv/T0;

% точность расчётов
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-40);
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-15);
% options = odeset('RelTol', 3e-14, 'AbsTol', 1e-11);
% options = odeset('RelTol', 3e-8, 'AbsTol', 1e-8);
if is_e==0
    options = odeset('RelTol', 1e-20, 'AbsTol', 1e-40);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-15);
end
xspan_m=[0  1e2
         0  1e2
         0  1e2    %5.4e-8 1e-6 % 3e-1
         0  1e2];
xspan=xspan_m(ind, :)/tau;
disp(xspan(2)*tau)
fileData = fopen('print_data_test.txt','w');
fileKrate = fopen('print_krate.txt','w');
% fprintf(fileData, '%% t     nO2(X)    nO2(a)  nO2(b)  nO(3P)  nO(1D)    O3  ne  T\n');
fprintf(fileData, '%% t     nO2(X)    nO2(a)  nO2(b)  nO2+    nO(3P)  nO(1D) nO+  nO- O3  O3exc   ne  T\n');
load LoKI-master\Code\Output\RNF\data4Maksim.mat
data4me=data4Maksim;           
pr_vec='%t';
for i=1:length(data4me.rate)
    pr_vec=[pr_vec '\t' data4me.rate(i).collDescription];
end
for i=1:length(data4me.rateExtra)
    pr_vec=[pr_vec '\t' data4me.rateExtra(i).collDescription];
end
fprintf(fileKrate, [pr_vec '\n']);
run_LoKI(T, n/n0, y(num_O2a), y(num_O2b), ne/n0, nO/n0, ...
                                                y(num_O1D), n0, EdN, 0);
                                        
disp('start')
% disp([y(num_T)])% y(num_E_e)])
tic
[X, Y]= ode15s(@(t,y) Rpart_Tiago_imitation(t, y, Prcl, Coll, tau, ...
    xspan, T0, n0, is_e, EdN, fileData, fileKrate),...
    xspan, y, options);             % рекомендуется использовать ode15s
X=X*tau;
Y(:, 1:num_e)=Y(:, 1:num_e)*n0;
Y(:, num_T)=Y(:, num_T)*T0;
data_Tiago_test=[X Y];

% run_LoKI(T, Y(end,num_O2)/n0, Y(end, num_O2a)/n0, Y(end, num_O2b)/n0, ...
%     ne/n0, Y(end, num_O)/n0, ...
%                                     Y(end, num_O1D)/n0, n0, EdN, 1e2);

toc
fclose('all');

%%
temp=data_Tiago_test;
e_i=O2.e_i(1,:);
num_vibr_levels=O2.num_vibr_levels(1);
% zero_p=sum(temp(1, 2:num_vibr_levels+1))*temp(1,end);
n_o2=sum(temp(:, 2:num_O2p+1), 2);
n_o=temp(:, num_O+1);
n_o1D=temp(:, num_O1D+1);
n_op=temp(:, num_Op+1);
n_om=temp(:, num_Om+1);
n_o3=0;
n_o3=temp(:, num_O3+1);
n_o3exc=temp(:, num_O3exc);
% all_p=n_o2.*temp(:,end);
rho=n_o2*O2.mass+(n_o+n_o1D+n_op+n_om)*O.mass+(n_o3+n_o3exc)*O3.mass;
disp(['max error rho ' num2str(max(abs((rho(1)-rho(end))/rho(1))))])
% tr vibr form rot
T_vec=temp(:, num_T+1);
% Tv = O2.e_i(1,2)./(k*log(temp(:, 2)./temp(:, 3)));
Energy=1.5*(n_o2+n_o+n_o1D+n_o3)*k.*T_vec+...
    sum(([O2.e_i(1,1:O2.num_vibr_levels(1)) ...
...    O2.e_i(2, 1:O2.num_vibr_levels(2)) ...
...    O2.e_i(3,1:O2.num_vibr_levels(3))...
    ]+O2.e_0).*...
    temp(:,2:1+num_O2b),2)+ n_o*O.form_e+(n_o2+n_o3)*k.*T_vec;%+...
%     1.5*ne*k*Tv;
disp(['max error Energy ' ...
    num2str(max(abs((Energy(1)-Energy(end))/Energy(1))))])
%% added processes:
%   O2(X)
% O2(X)  +O2(X)      ->2O(3P) +O2(X)    % diss
% O2(X)  +O2(a)      ->2O(3P) +O2(a)    % diss
% O2(X)  +O2(b)      ->2O(3P) +O2(b)    % diss
% O2(X)  +O(3P)      ->2O(3P) +O(3P)    % diss
% O2(a)  +O2(X)      ->2O(3P) +O2(X)    % diss
% 2O(3P) +O2(X)      ->O2(X)  +O2(X)    % rec
% 2O(3P) +O2(a)      ->O2(X)  +O2(a)    % rec
% 2O(3P) +O2(b)      ->O2(X)  +O2(b)    % rec
% 2O(3P) +O2(X)      ->O2(a)  +O2(X)    % rec
% 2O(3P) +O2(X)      ->O2(b)  +O2(X)    % rec
% 2O(3P) +O(3P)      ->O2(X)  +O(3P)    % rec
% O2(b)  +O(3P)      ->O2(X)  +O(3P)    % deex
% O2(a)  +O2(X)      ->O2(X)  +O2(X)    % deex
% O2(a)  +O2(a)      ->O2(X)  +O2(a)    % deex
% O2(a)  +O2(b)      ->O2(X)  +O2(b)    % deex
% O(1D)  +O2(X)      ->O(3P)  +O2(a)    % mix (deex, exci)
% O(1D)  +O2(X)      ->O(3P)  +O2(b)    % mix (deex, exci)
% O(1D)  +O2(X)      ->O(3P)  +O2(X)    % deex
% O2(a)  +O2(a)      ->O2(b)  +O2(X)    % exch
% O2(a)  +O(3P)      ->O2(X)  +O(3P)    % deex
% O(1D)  +O3(X)      ->2O2(X)           % exch (low)
% O(1D)  +O3(X)      ->O2(X)  +2O(3P)   % exch (low)
% O(3P)  +O2(X)+O(3P)->O3(X)  +O(3P)    % rec
% O(3P)  +O2(X)+O2(X)->O3(X)  +O2(X)    % rec
% O(3P)  +O2(X)+O3(X)->O3(X)  +O3(X)    % rec (low)
% O2(a)  +O3(X)      ->2O2(X) +O(3P)    % mix (diss, deex)
% O2(b)  +O3(X)      ->2O2(X) +O(3P)    % mix (diss, deex)
% O(3P)  +O3(X)      ->2O2(X)           % exch
% O(3P)  +O3(X)      ->O2(a)  +O2(X)    % exch
% O(3P)  +O3(X)      ->O2(b)  +O2(X)    % exch
% O(+)   +O3(X)      ->O2(+,X)+O2(X)    % exch
% O(+)   +O2(X)      ->O2(+,X)+O(3P)    % exch
% O(-)   +O(3P)      ->O2(X)  + e       % mix (rec, deion)
% O(-)   +O2(X)      ->O3(X)  + e       % mix (rec, deion)
% O(-)   +O2(b)      ->O(3P)  +O2(X)+e  % mix (deion, deexci)
% O2(+,X)+O(-)       ->O2(X)  +O(3P)    % deion
% O(3P)  +O2(X)+O2(X)->O3(exc)+O2(X)    % rec
% O3(exc)+O2(X)      ->O3(X)  +O2(X)    % deex
% O2(a)  +O3(exc)    ->2O2(X) +O(3P)    % mix (diss, deex)
% O(3P)  +O3(exc)    ->2O2(X)           % exch
% O2(X)  + e         ->2O(3P) + e       % elec, diss
% O2(X)  + e        <->O2(a)  + e       % elec, ex
% O2(X)  + e        <->O2(b)  + e       % elec, ex
% O2(X)  + e         ->O(3P)  +O(1D)+e  % elec, diss
% O3(X)  + e         ->O(3P)  +O2(X)+e  % elec, diss
% O2(X)  + e         ->O2(+,X)+ 2e      % elec, ion
% O2(X)  + e         ->O(3P)  + O(+)+2e % elec, mix (diss, ion)
% 	O2(a)
% O2(a) +O2(a)->2O(3P) +O2(a)    % diss
% O2(a) +O2(b)->2O(3P) +O2(b)    % diss
% O2(a) +O(3P)->2O(3P) +O(3P)    % diss 
% 2O(3P)+O2(a)->O2(a)  +O2(a)    % rec 
% 2O(3P)+O2(b)->O2(a)  +O2(b)    % rec
% 2O(3P)+O2(a)->O2(b)  +O2(a)    % rec
% O2(b) +O(3P)->O2(a)  +O(3P)    % deex
% O(+)  +O2(a)->O2(+,X)+O(3P)    % exch
% O(-)  +O2(a)->O3(X)  + e       % mix (rec, ion)
% O2(a) + e   ->2O(3P) + e       % elec, diss
% O2(a) + e  <->O2(b)  + e       % elec, ex
% O2(a) + e   ->O(3P)  +O(1D)+e  % elec, diss
% O2(a) + e   ->O2(+,X)+ 2e      % elec, ion
% O2(a) + e   ->O(3P)  +O(+) +2e % elec, mix (diss, ion)
% O2(X) + e   ->O(-)   +O(3P)    % elec, mix (diss, ion)
% O2(a) + e   ->O(-)   +O(3P)    % elec, mix (diss, ion)
%   O2(b)
% O2(b) +O2(b)->2O(3P)+O2(b)    % diss
% O2(b) +O(3P)->2O(3P)+O(3P)    % diss
% 2O(3P)+O2(b)->O2(b) +O2(b)    % rec
% O2(b) + e   ->2O(3P)+ e       % elec, diss
% O2(b) + e   ->O(3P) +O(1D)+e  % elec, diss
%   O(3P)
% O(1D)  +O(3P)->O(3P) + O(3P)   % deex
% O(+)   +O(-) ->2O(3P)          % deion
% O3(exc)+O(3P)->O3(X) + O(3P)   % deex
% O(3P)  + e  <->O(1D) +e        % elex, exci
% O2(+,X)+ e   ->2O(3P)          % elec, diss
% O2(+,X)+ e   ->O(3P) + O(1D)   % elec, diss
% O(3P)  + e   ->O(+)  + 2e      % elec, ion
% O(-)   + e   ->O(3P) + 2e      % elec, deion