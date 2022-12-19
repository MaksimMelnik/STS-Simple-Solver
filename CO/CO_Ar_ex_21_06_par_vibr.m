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

tic
for i_ini=[8]  % initial conditions test cases
                % 2 -- Fairbairn 8a,    3 -- Aliat's test case,
                % 4 -- Mick's fig 3a,   5 -- Appleton fig 3.
                % 6 -- mod Aliat low p, 7 -- mod Fairbairn 8a, high p
                % 8 -- Fairbairn 8f,    9 -- Fairbairn 8e
%  for i_exc=0%:1   % is electronic excitation on?
  for i_dis=1     % 1 -- Marrone-Treanor wo e exc., 2 -- MT with E exc.,
                        % 3 -- Aliat model, 4 -- Savelev model
   for i_U=3      % U parameter in MT and Aliat models
                        % 2 -- D/6k, 3 -- 3T, 4 -- Inf
%% variables 
% particles_data_ini;   % initialisation of particles and collisions data
load par_data.mat       % load of saved particle data
% is electronic excitaion on?
 ind_exc=1;      % COa and COA included (1), only CO(X) (0)
 if i_dis==1
  ind_exc=0;
 end
    V_K = 1.380649e-23;         % Boltzmann constant, J/K
    V_PI = 3.14159265358979323846;
%     checkX=1e-3;  % был такой параметр для вывода промежуточных значений

    % number of electronic states
if ind_exc==0
 CO.num_elex_levels=1;
end

    % ñîõðàíÿåì âñå ÷àñòèöû, ÷òîáû ïåðåäàòü ñêîïîì
Prcl.CO=CO; Prcl.C=C; Prcl.O=O; Prcl.C2=C2; Prcl.Ar=Ar; 
    
    % ñîõðàíÿåì âñå ñòîëêíîâåíèÿ, ÷òîáû ïåðåäàòü ñêîïîì
    Coll.CO_CO=Coll_CO_CO; Coll.CO_C=Coll_CO_C; Coll.CO_O=Coll_CO_O;
    Coll.CO_C2=Coll_CO_C2; Coll.CO_Ar=Coll_CO_Ar; Coll.C2=Coll_C2;
    Coll.CO_C__C2_O=Coll_CO_C__C2_O;
    
    Torr = 133.322368;
% ïåðåîïðåäåëÿåì âñ¸ ïîä Ôåéðáàðíà
        %   f       p0, Torr     v0, m/sec       T0, K      T2, K
init_c=[    
                % 1, âûñîêîòåìïåðàòóðíûé ñëó÷àé äëÿ ÷èñòîãî ÑÎ
            1       0.8         4746.548838     299.1113007914224 -1
                % 2, ñëó÷àé ñ àðãîíîì íà áàçå 8a Fairbairn. Äëÿ
                %   ñåðü¸çíûõ ìîäåëèðîâàíèé òåìïåðàòóðó âû÷èñëèòü
                %   ïåðåä ÓÂ
            0.01    10          2860            317.4471    7600
                % 3, ñëó÷àé ñòàòüè Àëèàòà
            1       500/Torr    5200            300         -1
                % 4, ïûòàåìñÿ ïîâòîðèòü òåñòêåéñ Ìèêà äëÿ ôèã3à
            0.001   5.7         2735            260         -1
                % 5, Appleton case for fig 3, T is approximate
            0.01    10          2090            300         -1
                % 6, ñëó÷àé ñòàòüè Àëèàòà, íî ïîìåíüøå ÷èñëîâàÿ
                %   ïëîòíîñòü
            1       0.04        5200            300         -1
                % 7, ñëó÷àé ñ àðãîíîì íà áàçå 8a Fairbairn, íî ñèëüíî
                %   áîëüøå äàâëåíèå
            0.01    100         2860            320         7600
                % 8, Fairbairn, случай 8f (T0 не точно)
            0.1     5.4         3260            300         9000
                % 9, Fairbairn, случай 8e (T0 не точно)
            0.001   10          3080            297         8800
        ];
%% easy access variables
num_COa=CO.num_vibr_levels(1)+1;
num_COA=num_COa+CO.num_vibr_levels(2);
num_C=num_COA+CO.num_vibr_levels(3);
num_O=num_C+1;
num_C2=num_O+1;
num_Ar=num_C2+1;
num_v=num_Ar+1;
num_T=num_v+1;
%%
% parfor ind_c=1:5
% for ind_c=1:1
  ind_c=i_ini;
 f=init_c(ind_c, 1);
 p0=init_c(ind_c, 2)*Torr;
 v0=init_c(ind_c, 3);

    T0=init_c(ind_c, 4);
    Tv0=T0;
    n0=p0/V_K/T0;
    NN = in_con_Ar([CO.mass, v0, T0, Ar.mass, f]);
    n1 = NN(1);     % áåçðàçìåðíûå
    T1 = NN(2);     % áåçðàçìåðíûå
    v1 = NN(3);     % áåçðàçìåðíûå
    disp([num2str(ind_c) ': T1=' num2str(T1*T0) ', T0=' num2str(T0) ...
            ', v1=' num2str(v1*v0) ', n1=' num2str(n1*n0, '%1.3e')])
%         pause
    ys=zeros(num_T,1);
    n=density_f_exc(Tv0, n1*f, CO);
    ys(1:num_COA+CO.num_vibr_levels(3)-1)=n;
    ys(num_C)=0;%2e-1;                                  % C
    ys(num_O)=0;%2e-1;                                  % O
    ys(num_C2)=0;                                       % C2
    ys(num_Ar)=n1*(1-f);                                % Ar
    ys(num_v)=v1;                                       % v
    ys(num_T)=T1;                                       % T
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-30, ...
%                                     'NonNegative', 1:num_T); %#ok<NASGU>
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-18, ...
%                                     'NonNegative', 1:num_T); %#ok<NASGU>
%     options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-13, ...
%                                     'NonNegative', 1:num_T); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-11, ...
%                                     'NonNegative', 1:num_T); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-8, ...
%                                     'NonNegative', 1:num_T); %#ok<NASGU>
%     options_s = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, ...
%                                                 'NonNegative', 1:num_T);
    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, ...
                                    'NonNegative', 1:num_T); %%#ok<NASGU>
%     options_s = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, ...
%                                                 'NonNegative', 1:num_T);
    sigma0 = V_PI*Coll_CO_CO.coll_diameter^2;
    Delta = 1 / sqrt(2) / n0 / sigma0;
 xspan=[0 1e0]/Delta; % 1e5, 4e1 for tc7 test
    disp([num2str(ind_c) ': xspan=' num2str(xspan(2),'%.2e')])
    ind_Arr=1;      % Arrenius law constants: 1 - Park, 2 - Ibraguimova,
                        % 3 - McKenzie, 4 -- Fairbairn, 5 -- Appleton
                        % 6 -- Appleton with changed diss energy, 
                        % 7 -- Mick
	struct_Arr(1).text='Park';      struct_Arr(2).text='Ibraguimova';
    struct_Arr(3).text='McKenzie';  struct_Arr(4).text='Fairbairn';
    struct_Arr(5).text='Appleton (Arrhenius)';  
    struct_Arr(6).text='Appleton changed diss E';
    struct_Arr(7).text='Mick';
	ind_U=i_U;
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
 setup.f=f;
 switch setup.model_VT
     case 'FHO'
         disp('VT rates by old FHO code.')
     case 'SSH'
         disp('VT rates by SSH theory.')
 end
 if setup.C2
     disp('C2 is included.')
 end
 [X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, n0, T0, v0, Delta, Prcl, Coll, ...
               ind_Arr, ind_U, i_dis, ind_exc, setup),...
                                                xspan, ys, options_s);
 X=X*Delta;
 Y(:,num_T)=Y(:,num_T)*T0;
 Y(:,num_v)=Y(:,num_v)*v0;
 Y(:,1:num_Ar)=Y(:,1:num_Ar)*n0;
 T=Y(:, num_T);
 Tv = CO.ev_i{1}(2)./(V_K*log(Y(:,1)./Y(:,2)));
 time_ms=X./v0*1e6;
 out=[X, Y, Tv, time_ms];
 temp=out;
 res(i_ini, i_dis, i_U).temp=temp; %#ok<AGROW>
    
 nCO=sum(temp(:,2:1+CO.num_vibr_levels(1)),2);                   % n_CO
 rho=(nCO+sum(Y(:,num_COa:num_COa+CO.num_vibr_levels(2)-1),2)... % rho
    +sum(Y(:,num_COA:num_COA+CO.num_vibr_levels(3)-1),2))*CO.mass...
    +Y(:, num_C)*C.mass + Y(:,num_O)*O.mass + Y(:, num_C2)*C2.mass...
                                                    +Y(:,num_Ar)*Ar.mass;
 pres=sum(Y(:,1:num_Ar),2)*V_K.*T;                           % p
 rhov=rho.*Y(:, num_v);                                      % rho*v
 disp([num2str(ind_c) ': rho*v max error ' ...
                            num2str(max(abs((rhov-rhov(1))/rhov(1))))])
 rhov2p=rho.*Y(:, num_v).^2+pres;                            % rho*v^2+p
 disp([num2str(ind_c) ': rho*v^2+p max error '...
                        num2str(max(abs((rhov2p-rhov2p(1))/rhov2p(1))))])
 En=2.5*(sum(Y(:,1:num_COA+CO.num_vibr_levels(3)-1), 2)...    % E
                                                +Y(:,num_C2))*V_K.*T...
    +1.5*(Y(:,num_C)+Y(:,num_O)+Y(:,num_Ar))*V_K.*T...
    +Y(:,1:CO.num_vibr_levels(1))*(CO.ev_i{1}+CO.ev_0(1))'...
    +Y(:,num_COa:num_COa+CO.num_vibr_levels(2)-1)...
                                *(CO.ev_i{2}+CO.ev_0(2)+CO.e_E(2))'...
	+Y(:,num_COA:num_COA+CO.num_vibr_levels(3)-1)...
                                *(CO.ev_i{3}+CO.ev_0(3)+CO.e_E(3))'...
    +sum(Y(:,1:num_COA+CO.num_vibr_levels(3)-1),2)*CO.form_e...
    +Y(:,num_C)*C.form_e+Y(:,num_O)*O.form_e;
 Ep=(En+pres)./rho+0.5*Y(:, num_v).^2;                  % E conservation
 disp([num2str(ind_c) ': E cons max error ' ...
                                    num2str(max(abs((Ep-Ep(1))/Ep(1))))])
   end
  end
%  end
end
 figure
 semilogx(temp(:,1), T, temp(:,1), Tv)
out=res;
toc
end