function out=CO_Ar_ex_21_06_par_vibr
% Программа для расчёта смесей CO и Ar с учётом CO* и C2 с целью выявить
% задержку диссоциации. На базе кода для параллельных расчётов смеси CO/Ar
% CO_Ar_mixure_20_04_par.m, модификация кода CO_Ar_ex_20_10_par_vibr с 
% учётом диссоциации по Алиату и колебательных уровней CO*, + диссоциация
% по Савельеву и прочие дополнения для Алушты
% 23.06.2021

%% warnings
% warning('VE processes off.')
% warning('xspan for test')
% warning('kVE down in kVE and RVE: CO(a,A)->CO')
% warning("kVE UP in kVE and RVE, but Aliat's transitions only.")
% warning("kVE DOWN in kVE and RVE, but Aliat's transitions only.")
% warning('Calculation accuracy is under a test.')
% warning('CO+C->C2+O is under the test')

    % constants
N_a=6.02214076e23;          % Avogadro constant
k=1.380649e-23;             % Boltzmann constant, J/K
Torr = 133.322368;          % how much Pa in Torr
addpath('../CO/')

tic
for i_ini=[8]  % initial conditions test cases
                % 2 -- Fairbairn 8a,    3 -- Aliat's test case,
                % 4 -- Mick's fig 3a,   5 -- Appleton fig 3.
                % 6 -- mod Aliat low p, 7 -- mod Fairbairn 8a, high p
                % 8 -- Fairbairn 8f,    9 -- Fairbairn 8e
%  for i_exc=0%:1   % is electronic excitation on?
  for i_dis=3     % 1 -- Marrone-Treanor wo e exc., 2 -- MT with E exc.,
                        % 3 -- Aliat model, 4 -- Savelev model
   for ind_U=3      % U parameter in MT and Aliat models
                        % 2 -- D/6k, 3 -- 3T, 4 -- Inf
%% variables 
% particles_data_ini;   % initialisation of particles and collisions data
load par_data.mat       % load of saved particle data
% is electronic excitaion on?
 ind_exc=1;      % COa and COA included (1), only CO(X) (0)
 if i_dis==1
  ind_exc=0;
 end
%     checkX=1e-3;  % был такой параметр для вывода промежуточных значений

    % number of electronic states
if ind_exc==0
 CO.num_elex_levels=1;
end
    
        %   f       p0, Torr     v0, m/sec       T0, K      T2, K
init_c=[    
                % 1, высокотемпературный случай для чистого СО
            1       0.8         4746.548838     299.1113007914224 -1
                % 2, случай с аргоном на базе 8a Fairbairn. Для
                %   серьёзных моделирований температуру вычислить
                %   перед УВ
            0.01    10          2860            317.4471    7600
                % 3, случай статьи Алиата
            1       500/Torr    5200            300         -1
                % 4, пытаемся повторить тесткейс Мика для фиг3а
            0.001   5.7         2735            260         -1
                % 5, Appleton case for fig 3, T is approximate
            0.01    10          2090            300         -1
                % 6, случай статьи Алиата, но поменьше числовая
                %   плотность
            1       0.04        5200            300         -1
                % 7, случай с аргоном на базе 8a Fairbairn, но сильно
                %   больше давление
            0.01    100         2860            320         7600
                % 8, Fairbairn, случай 8f (T0 не точно)
            0.1     5.4         3260            300         9000
                % 9, Fairbairn, случай 8e (T0 не точно)
            0.001   10          3080            297         8800
        ];
%%
% parfor ind_c=1:5
% for ind_c=1:1
 f=init_c(i_ini, 1);
 p0=init_c(i_ini, 2)*Torr;
 v0=init_c(i_ini, 3);
 T0=init_c(i_ini, 4);
    Tv0=T0;
    n0=p0/k/T0;
    NN = in_con_Ar([CO.mass, v0, T0, Ar.mass, f]);
    n1 = NN(1);     % dimensionless
    T1 = NN(2);     % dimensionless
    v1 = NN(3);     % dimensionless
 disp([num2str(i_ini) ': T1=' num2str(T1*T0) ', T0=' num2str(T0) ...
                ', v1=' num2str(v1*v0) ', n1=' num2str(n1*n0, '%1.3e')])
    sigma0 = pi*Coll_CO_CO.coll_diameter^2;
    Delta = 1 / sqrt(2) / n0 / sigma0;
 xspan=[0 1e0]/Delta; % 1e5, 4e1 for tc7 test
    disp([num2str(i_ini) ': xspan=' num2str(xspan(2),'%.2e')])
    ind_Arr=1;      % Arrenius law constants: 1 - Park, 2 - Ibraguimova,
                        % 3 - McKenzie, 4 -- Fairbairn, 5 -- Appleton
                        % 6 -- Appleton with changed diss energy, 
                        % 7 -- Mick
	struct_Arr(1).text='Park';      struct_Arr(2).text='Ibraguimova';
    struct_Arr(3).text='McKenzie';  struct_Arr(4).text='Fairbairn';
    struct_Arr(5).text='Appleton (Arrhenius)';  
    struct_Arr(6).text='Appleton changed diss E';
    struct_Arr(7).text='Mick';
	struct_U(1).text='Savelev?'; struct_U(2).text='D/6k';
    struct_U(3).text='3T'; struct_U(4).text='Inf';
 struct_Aliat(1).text='Marrone-Treanor'; 
 struct_Aliat(2).text='Marrone-Treanor'; struct_Aliat(3).text='Aliat';
 struct_Aliat(4).text='Savelev for CO+CO, Aliat for other collisions';
 disp(['Arrhenius dissociation constants are provided by ' ...
                                           struct_Arr(ind_Arr).text '.'])
 disp(['U parameter in non-eq dissociation model is ' ...
                                               struct_U(ind_U).text '.'])
 disp(['Dissociation model by ' struct_Aliat(i_dis).text '.'])
 disp(['CO(X)' ind_exc*', CO(a) and CO(A)' ' states are included.'])
 setup.C2=1;        % is C2 included?
 setup.model_VT='FHO';
 setup.model_VT='SSH';
 switch setup.model_VT
     case 'FHO'
         disp('VT rates by old FHO code.')
     case 'SSH'
         disp('VT rates by SSH theory.')
 end
 if setup.C2
     disp('C2 is included.')
 end
 
switch ind_Arr
    case 1
        Diss.Arrhenius='Park';
        keys={'CO', 'C', 'O', 'Ar', 'C2'};
        val_A={2.3e20/N_a*1e-6, 3.4e20/N_a*1e-6, 3.4e20/N_a*1e-6, ...
                                    2.3e19/N_a*1e-6, 2.3e20/N_a*1e-6};
        val_n={-1, -1, -1, -1, -1};
        diss_Arrhenius_A_CO=containers.Map(keys, val_A);
        diss_Arrhenius_n_CO=containers.Map(keys, val_n);
        val_A={3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6, ...
                                    3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6};
        val_n={0, 0, 0, 0, 0};
        diss_Arrhenius_A_C2=containers.Map(keys, val_A);
        diss_Arrhenius_n_C2=containers.Map(keys, val_n);
    case 2
        Diss.Arrhenius='Ibraguimova';
    case 3
        Diss.Arrhenius='McKenzie';
    case 4
        Diss.Arrhenius='Fairbairn';
end
mix.num=0;
CO.diss_Arrhenius_A=diss_Arrhenius_A_CO;
CO.diss_Arrhenius_n=diss_Arrhenius_n_CO;
Ps={mix, CO, C, O};
if setup.C2
    C2.diss_Arrhenius_A=diss_Arrhenius_A_C2;
    C2.diss_Arrhenius_n=diss_Arrhenius_n_C2;
    Ps{length(Ps)+1}=C2;
end
if f<1
    Ps{length(Ps)+1}=Ar;
end
num=0;
index{1}=0;
for ind=2:length(Ps)
    num_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_states-1;
end
Diss.rec=true;              % Is recombination included?
switch i_dis
    case 1
        Diss.NEmodel='MT';
    case 2
        Diss.NEmodel='MT';
    case 3
        Diss.NEmodel='Aliat';
    case 4
        Diss.NEmodel='Savelev21';
end
switch ind_U
    case 2
        Diss.U='D/6k';
    case 3
        Diss.U='3T';
    case 4
        Diss.U='inf';
end
Reacs_keys={'Diss', 'VT'};  % Exch?
reacs_val={Diss, setup.model_VT};
if ind_exc
    Reacs_keys=[Reacs_keys, 'VE'];
    reacs_val=[reacs_val, 'VE'];
end
kinetics.Ps=Ps(2:end);
kinetics.num_Ps=length(kinetics.Ps);
kinetics.num_eq=num;
kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
kinetics.index=index(2:end);
kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;

 ys2=zeros(kinetics.num_eq+2, 1);
 n=density_f_exc(Tv0, n1*f, CO);
 ys2(1:length(n))=n;
 if f<1
  ys2(kinetics.index{end})=n1*(1-f);
 end
 ys2(end-1)=v1;
 ys2(end)=T1;
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-30, ...
%                         'NonNegative', 1:kinetics.num_eq+2); %#ok<NASGU>
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-18, ...
%                         'NonNegative', 1:kinetics.num_eq+2); %#ok<NASGU>
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-13, ...
%                         'NonNegative', 1:kinetics.num_eq+2); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-11, ...
%                         'NonNegative', 1:kinetics.num_eq+2); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-8, ...
%                         'NonNegative', 1:kinetics.num_eq+2); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, ...
%                                     'NonNegative', 1:kinetics.num_eq+2);
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, ...
                        'NonNegative', 1:kinetics.num_eq+2); %%#ok<NASGU>
%     options_s = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, ...
%                                     'NonNegative', 1:kinetics.num_eq+2);
 [X2, Y2]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics),...
                                                xspan, ys2, options_s);
                                
 X2=X2*Delta;
 Y2(:, end-1) = Y2(:, end-1)*v0;
 Y2(:, end) =   Y2(:, end)*T0;
 Y2(:, 1:end-2)=Y2(:, 1:end-2)*n0;
 T=Y2(:, end);
 Tv = CO.ev_i{1}(2)./(k*log(Y2(:,1)./Y2(:,2)));
 time_ms=X2./v0*1e6;
 res(i_ini, i_dis, ind_U).temp=[X2, Y2, Tv, time_ms]; %#ok<AGROW>
 
 rhov0=n0*(f*CO.mass+(1-f)*Ar.mass) * v0;                    % rho0*v0
 rhov2p0=n0*(f*CO.mass+(1-f)*Ar.mass) * v0^2 + p0;
 e_i=[];
 for ind_e=1:CO.num_elex_levels
  e_i=[e_i, CO.ev_i{ind_e}+CO.ev_0(ind_e)+CO.e_E(ind_e)];
 end
 En0=n0*e_i*n/n1+n0*f*k*T0 + 1.5*n0*k*T0+n0*f*CO.form_e;
 Ep0=(En0+p0)/(n0*(f*CO.mass+(1-f)*Ar.mass))+0.5*v0^2;
 disp([num2str(i_ini) ': conservation laws check'])
 check_CL_SW([rhov0 rhov2p0 Ep0], Y2, kinetics);
   end
  end
%  end
end
 figure
 semilogx(X2, T, X2, Tv)
out=res;
toc
rmpath('../CO/')
end