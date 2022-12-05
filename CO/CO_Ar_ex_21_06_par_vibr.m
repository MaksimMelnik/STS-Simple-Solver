function out=CO_Ar_ex_21_06_par_vibr
% Ïðîãðàììà äëÿ ðàñ÷¸òà ñìåñåé CO è Ar ñ ó÷¸òîì CO* è C2 ñ öåëüþ âûÿâèòü
% çàäåðæêó äèññîöèàöèè. Íà áàçå êîäà äëÿ ïàðàëëåëüíûõ ðàñ÷¸òîâ ñìåñè CO/Ar
% CO_Ar_mixure_20_04_par.m, ìîäèôèêàöèÿ êîäà CO_Ar_ex_20_10_par_vibr ñ 
% ó÷¸òîì äèññîöèàöèè ïî Àëèàòó è êîëåáàòåëüíûõ óðîâíåé CO*, + äèññîöèàöèÿ
% ïî Ñàâåëüåâó è ïðî÷èå äîïîëíåíèÿ äëÿ Àëóøòû
% 23.06.2021

%% warnings
% warning('Only ground electronic CO state in dissociation.')
% warning('VE processes off.')
% warning('Detailed balance is old and uncorrect. For C and O s_e = 1.')
% warning('Old FHO code not allowed for CO-Ar.')
% warning('dE in diss detailed balance and form_e probably uncorrect.')
% warning('FHO under the test')
% warning('xspan for test')
% warning('Dissociation off.')
% warning('Ar: dissociation in collision with Ar is off.')
% warning('kVE down in kVE and RVE: CO(a,A)->CO')
% warning('kVE UP in kVE and RVE, but Aliats transitions only.')
% warning('kVE DOWN in kVE and RVE, but Aliats transitions only.')
% warning('w/o exc states modification in progress')
% warning('Recombination in Aliat diss function is off.')
% warning('Recombination in Savelev diss function is off.')
% warning('Calculation accuracy is under a test.')
% warning('Changed diss energy (Appleton)')
% warning('Changed Arrenius law parameters')
% warning('C2 included')
% warning('kd is constant')
% warning('params for the laptop heat test')
% warning('CO+C->C2+O is under the test')

tic
for i_ini=[3]  % initial conditions test cases
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
% is electronic excitaion on?
 ind_exc=1;      % COa and COA included (1), only CO(X) (0)
 if i_dis==1
  ind_exc=0;
 end
%     V_H = 6.626070041e-34;      % Plank constant, J*sec
    V_K = 1.380649e-23;         % Boltzmann constant, J/K
    V_PI = 3.14159265358979323846; 
%     V_C = 299792458; % speed of light
%     V_HBAR = 1.054571800e-34; V_R = 8.3144621;
%     checkX=1e-3;  % áûë òàêîé ïàðàìåòð äëÿ âûâîäà ïðîìåæóòî÷íûõ çíà÷åíèé

    particles_data_ini;   % initialisation of particles and collisions 
                                    %   data
    load par_data.mat       % load of saved particle data

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
                % 8, Fairbairn, ñëó÷àé 8f (T0 íå òî÷íî)
            0.1     5.4         3260            300         9000
                % 8, Fairbairn, ñëó÷àé 8e (T0 íå òî÷íî)
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
 disp('VT probabilities obtained using probably old FHO code.')
 disp(['CO(X)' ind_exc*', CO(a) and CO(A)' ' states are included.'])
 setup.C2=1;        % is C2 included?
 [X, Y]=ode15s(@(t, y) Rpart_s(t, y, n0, T0, v0, Delta, Prcl, Coll, ...
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

    % äëÿ ÓÂ
function out=Rpart_s(t, y, n0, T0, v0, Delta, Prcl, Coll, ind_Arr,...
                                ind_U, i_dis, ind_exc, setup) %#ok<INUSL>
ind_Aliat=0;
if i_dis>2
    ind_Aliat=1;
end
% disp(t)
num_vibr_levels=Prcl.CO.num_vibr_levels(1);
num_COa=num_vibr_levels+1;
num_COA=num_COa+Prcl.CO.num_vibr_levels(2);
num_C=num_COA+Prcl.CO.num_vibr_levels(3);
num_O=num_C+1;
num_C2=num_O+1;
num_Ar=num_C2+1;
num_v=num_Ar+1;
num_T=num_v+1;
k=1.380648528e-23;
e_i =(Prcl.CO.ev_i{1}+Prcl.CO.ev_0(1))/k/T0; 
e_ia=(Prcl.CO.ev_i{2}+Prcl.CO.ev_0(2))/k/T0;
e_iA=(Prcl.CO.ev_i{3}+Prcl.CO.ev_0(3))/k/T0;
                                    % levels_e(num_vibr_levels)'/k/T0;
ni_b = y(1:num_vibr_levels);
nCOa_b = y(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
nCOA_b = y(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
nac_b = y(num_C);
nao_b = y(num_O);
nC2_b = y(num_C2);
naAr_b = y(num_Ar);
nm_b = sum(ni_b);
v_b = y(num_v);
T_b = y(num_T);

    % âûäåëÿåì ïàìÿòü
R=zeros(num_T,1);
COnl0=sum(Prcl.CO.num_vibr_levels);
R_d_CO_CO=zeros(COnl0,1);   R_d_CO_C=zeros(COnl0,1);        %#ok<PREALL>
R_d_CO_O=zeros(COnl0,1);    R_d_CO_C2=zeros(COnl0,1);       %#ok<PREALL>
R_d_CO_Ar=zeros(COnl0,1);                                   
R_d_COa_CO=0; R_d_COa_COa=0; R_d_COa_COA=0; R_d_COa_C=0;    %#ok<NASGU>
R_d_COa_O=0; R_d_COa_C2=0; R_d_COa_Ar=0;                    %#ok<NASGU>
R_d_COA_CO=0; R_d_COA_COa=0; R_d_COA_COA=0; R_d_COA_C=0;    %#ok<NASGU>
R_d_COA_O=0; R_d_COA_C2=0; R_d_COA_Ar=0;                    %#ok<NASGU>
R_d_C2=0;                                                   %#ok<NASGU>
R_exch_CO_C_C2_O=zeros(COnl0,1);
R_VT_CO_CO=zeros(num_vibr_levels,1);                        %#ok<PREALL>
R_VT_CO_C=zeros(num_vibr_levels,1);                         %#ok<PREALL>
R_VT_CO_O=zeros(num_vibr_levels,1);                         %#ok<PREALL>
R_VT_CO_C2=zeros(num_vibr_levels,1);
R_VT_CO_Ar=zeros(num_vibr_levels,1);                        
R_VT_data_COa=0;    R_VT_data_COA=0;                        
R_VE_CO_toCOa=zeros(num_vibr_levels,...
                            Prcl.CO.num_vibr_levels(2));    
R_VE_CO_toCOA=zeros(num_vibr_levels,1);                     
R_VE_dis=0;                                                 %#ok<NASGU>
R_d_CO_C2_A=0;                                              %#ok<NASGU>
    % Ïåðå÷åíü ó÷ò¸ííûõ ïðîöåññîâ. 
    % ×òîáû îòêëþ÷èòü, ïðîñòî çàêîììåíòèðîâàòü.
    % CO dissociation
if i_dis<4
    R_d_CO_CO=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], ...
            nac_b, nao_b, nm_b+sum(nCOa_b)+sum(nCOA_b), Coll.CO_CO, ...
                                 T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
else
 ind_exc_Savelev=3; % 1 - òîëüêî äèññîöèèðóþùèé àòîì âîçáóæä¸í
                    % 2 - åù¸ + ïàðòí¸ð, 3 - åù¸ + ïðîäóêò
    R_d_CO_CO_Savelev=R_diss_Savelev(Prcl.CO, [ni_b; nCOa_b; nCOA_b], ...
        nac_b, nao_b, nm_b+sum(nCOa_b)+sum(nCOA_b), Coll.CO_CO, ...
                            T_b*T0, n0, ind_Arr, ind_U, ind_exc_Savelev);
    R_d_CO_CO=sum(R_d_CO_CO_Savelev, [2, 3]);
%     R_VE_dis=-reshape(sum(R_d_CO_CO_Savelev, [1, 2]), [], 1, 1)...
%                                        +sum(R_d_CO_CO_Savelev, [1, 3])';
% for ind_exc_Savelev=2
% R_VE_dis=-sum(R_d_CO_CO_Savelev, [1, 3])';
% R_VE_dis(1)=R_VE_dis(1)+sum(R_d_CO_CO_Savelev, [1, 2]);
%     òàê, íåò, ind_Aliat = 0 èëè 1. Íóæíî ïåðåäàâàòü i_dis âñ¸ òàêè
end
R_d_CO_C= R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b, ...
        nao_b, nac_b, Coll.CO_C, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
R_d_CO_O= R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b, ...
        nao_b, nao_b, Coll.CO_O, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
if naAr_b>0
 R_d_CO_Ar=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b,...
      nao_b, naAr_b, Coll.CO_Ar, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
end
if setup.C2
 R_d_CO_C2=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b,...
      nao_b, nC2_b, Coll.CO_C2, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
end
R_d_data=R_d_CO_CO+R_d_CO_C+R_d_CO_O+R_d_CO_C2+R_d_CO_Ar;
    % COa dissociation
R_d_COa=R_d_data(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
    % COA dissociation
% R_d_COA=0;
R_d_COA=R_d_data(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
    % C2 dissociation
% R_d_C2=R_diss(Prcl.C2, nC2_b, nac_b, nac_b, ...   fix s_e in Kdr
%     nm_b+nCOa_b+nCOA_b+nac_b+nao_b+nC2_b+naAr_b,...
%     Coll.C2, T_b*T0, n0, ind_Arr, ind_U);
if setup.C2
 R_d_C2=R_diss_Aliat_onoff(Prcl.C2, nC2_b, nac_b, nac_b, ...
      nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, ...
                        Coll.C2, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
 R_exch_CO_C_C2_O=R_exch_CO_C__C2_O(Prcl.C2, Prcl.CO, ...
     [ni_b; nCOa_b; nCOA_b], nac_b, nao_b, nC2_b, Coll.CO_C__C2_O, ...
                                                       T_b*T0, ind_Arr);
end
    % CO VT
% Try to use integration for FHO kVT calculation. Differs from the kVT
%   rate coefficients in Aliat 2003 paper.
% R_VT_CO_CO =R_VT_exc(Prcl.CO, ni_b, Prcl.CO, ...
%                             nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 1);
% R_VT_CO_C =R_VT_exc(Prcl.CO, ni_b, Prcl.C, nac_b, T_b*T0, 1);
% R_VT_CO_O =R_VT_exc(Prcl.CO, ni_b, Prcl.O, nao_b, T_b*T0, 1);
% R_VT_CO_C2=R_VT_exc(Prcl.CO, ni_b, Prcl.C2, nC2_b, T_b*T0, 1);
% R_VT_CO_Ar=R_VT_exc(Prcl.CO, ni_b, Prcl.Ar, naAr_b, T_b*T0, 1);

% Old version of FHO kVT using solver for non-linear equations. Code
% provided by Olga Kunova.
R_VT_CO_CO=R_VT_old(Prcl.CO, ni_b, Prcl.CO, ...
                              nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 1);
R_VT_CO_C =R_VT_old(Prcl.CO, ni_b, Prcl.C,  nac_b,  T_b*T0, 1);
R_VT_CO_O =R_VT_old(Prcl.CO, ni_b, Prcl.O,  nao_b,  T_b*T0, 1);
if naAr_b>0
 R_VT_CO_Ar=R_VT_old(Prcl.CO, ni_b, Prcl.Ar, naAr_b, T_b*T0, 1);
end
if setup.C2
 R_VT_CO_C2=R_VT_old(Prcl.CO, ni_b, Prcl.C2, nC2_b, T_b*T0, 1);
end

R_VT_data=R_VT_CO_CO+R_VT_CO_C+R_VT_CO_O+R_VT_CO_C2+R_VT_CO_Ar;
% R_VT_data_COa=R_VT_exc(Prcl.CO, nCOa_b, Prcl.CO, ...        CO, COa, COA
%                             nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 2)...
%            +R_VT_exc(Prcl.CO, nCOa_b, Prcl.C, nac_b, T_b*T0, 2)...   C
%            +R_VT_exc(Prcl.CO, nCOa_b, Prcl.O, nao_b, T_b*T0, 2)...   O
%            ...+R_VT_exc(Prcl.CO, nCOa_b, Prcl.C2, nC2_b, T_b*T0, 2)...C2
%            +R_VT_exc(Prcl.CO, nCOa_b, Prcl.Ar, naAr_b, T_b*T0, 2)... Ar
%            ;
% pause
if ind_exc
 R_VT_data_COa=R_VT_old(Prcl.CO, nCOa_b, Prcl.CO, ...         CO, COa, COA
                            nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 2)...
           +R_VT_old(Prcl.CO, nCOa_b, Prcl.C, nac_b, T_b*T0, 2)...   C
           +R_VT_old(Prcl.CO, nCOa_b, Prcl.O, nao_b, T_b*T0, 2)...   O
           ;
 if naAr_b>0
  R_VT_data_COa=R_VT_data_COa...
           +R_VT_old(Prcl.CO, nCOa_b, Prcl.Ar, naAr_b, T_b*T0, 2);
 end
 if setup.C2
  R_VT_data_COa=R_VT_data_COa+...
                    R_VT_old(Prcl.CO, nCOa_b, Prcl.C2, nC2_b, T_b*T0, 2);
 end
% R_VT_data_COA=R_VT_exc(Prcl.CO, nCOA_b, Prcl.CO, ...        CO, COa, COA
%                             nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 3)...
%            +R_VT_exc(Prcl.CO, nCOA_b, Prcl.C, nac_b, T_b*T0, 3)...   C
%            +R_VT_exc(Prcl.CO, nCOA_b, Prcl.O, nao_b, T_b*T0, 3)...   O
%            ...+R_VT_exc(Prcl.CO, nCOA_b, Prcl.C2, nC2_b, T_b*T0, 3)...C2
%            +R_VT_exc(Prcl.CO, nCOA_b, Prcl.Ar, naAr_b, T_b*T0, 3)... Ar
%            ;
 R_VT_data_COA=R_VT_old(Prcl.CO, nCOA_b, Prcl.CO, ...         CO, COa, COA
                            nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 3)...
           +R_VT_old(Prcl.CO, nCOA_b, Prcl.C, nac_b, T_b*T0, 3)...   C
           +R_VT_old(Prcl.CO, nCOA_b, Prcl.O, nao_b, T_b*T0, 3)...   O
           ;
 if naAr_b>0
  R_VT_data_COA=R_VT_data_COA...
           +R_VT_old(Prcl.CO, nCOA_b, Prcl.Ar, naAr_b, T_b*T0, 3);
 end
 if setup.C2
  R_VT_data_COA=R_VT_data_COA+...
                    R_VT_old(Prcl.CO, nCOA_b, Prcl.C2, nC2_b, T_b*T0, 3);
 end

    % CO VE
 R_VE_CO_toCOa=R_VE_m(Prcl.CO, ni_b, nCOa_b, ...
    nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, T_b*T0, 2);
 R_VE_CO_toCOA=R_VE_m(Prcl.CO, ni_b, nCOA_b, ...
    nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, T_b*T0, 3);
end

R(1:num_vibr_levels)=R_VT_data+R_d_data(1:num_vibr_levels)+...
    sum(R_VE_CO_toCOa,2)+sum(R_VE_CO_toCOA,2);
R(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1)=R_d_COa+R_VT_data_COa...
                                                -sum(R_VE_CO_toCOa,1)';
R(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1)=R_d_COA+R_VT_data_COA...
                                                -sum(R_VE_CO_toCOA,1)';
R(num_C)=-sum(R_d_data)-2*R_d_C2+sum(R_exch_CO_C_C2_O);
R(num_O)=-sum(R_d_data)-sum(R_exch_CO_C_C2_O);
R(num_C2)=R_d_C2-sum(R_exch_CO_C_C2_O);
R(1:num_COA+Prcl.CO.num_vibr_levels(3)-1)=...
    R(1:num_COA+Prcl.CO.num_vibr_levels(3)-1)+R_VE_dis+R_exch_CO_C_C2_O;
% disp(sum(R(1:num_C-1)*Prcl.CO.mass)+R(num_C)*Prcl.C.mass...
%     +R(num_O)*Prcl.O.mass+R(num_Ar)*Prcl.Ar.mass)
% disp(R(68)*Prcl.CO.mass)
% CO=Prcl.CO; C=Prcl.C; O=Prcl.O;
% disp(sum(R_d_CO_CO)*CO.mass- sum(R_d_CO_CO)*(C.mass+O.mass))
% disp(sum(R_d_CO_CO+R_d_CO_C+R_d_CO_O)*CO.mass...
%     - sum(R_d_CO_CO+R_d_CO_C+R_d_CO_O)*(C.mass+O.mass))
% disp(sum(R_VT_CO_CO))
% pause

% if t>0
%     disp(R(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1)')
%     disp('VE')
%     disp(-sum(R_VE_CO_toCOA))
%     pause
% end

R=R*n0*Delta/v0;

    % íåðàçðûâíîñòü
M=zeros(num_T);
for i=1:num_Ar
    M(i,i)=v_b;
end
% M(1:num_vibr_levels,num_v)=ni_b';
M(1:num_Ar, num_v)=y(1:num_Ar)';
    % èìïóëüñ
M(num_v,1:num_Ar)=T_b;
M(num_v,num_v)=(Prcl.CO.mass*(nm_b+sum(nCOa_b)+sum(nCOA_b))+...
    Prcl.C.mass*nac_b+Prcl.O.mass*nao_b+...
    Prcl.C2.mass*nC2_b+Prcl.Ar.mass*naAr_b)*v_b*v0^2/k/T0;
M(num_v,num_T)=nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b;
    % ýíåðãèÿ
M(num_T,1:num_vibr_levels) = 2.5*T_b+e_i+Prcl.CO.form_e/k/T0;
M(num_T,num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1)=...
                    2.5*T_b+e_ia+Prcl.CO.form_e/k/T0+Prcl.CO.e_E(2)/k/T0;
M(num_T,num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1)=...
                    2.5*T_b+e_iA+Prcl.CO.form_e/k/T0+Prcl.CO.e_E(3)/k/T0;
M(num_T,num_C)=1.5*T_b+Prcl.C.form_e/k/T0;
M(num_T,num_O)=1.5*T_b+Prcl.O.form_e/k/T0;
M(num_T,num_C2)=2.5*T_b;
M(num_T,num_Ar)=1.5*T_b;
M(num_T,num_v)=1/v_b*(3.5*(nm_b+sum(nCOa_b)+sum(nCOA_b)+nC2_b)*T_b...
    +2.5*(nac_b+nao_b+naAr_b)*T_b+e_i*ni_b+e_ia*nCOa_b+e_iA*nCOA_b...
    +(nm_b+sum(nCOa_b)+sum(nCOA_b))*Prcl.CO.form_e/k/T0...
    +nac_b*Prcl.C.form_e/k/T0+nao_b*Prcl.O.form_e/k/T0...
    +sum(nCOa_b)*Prcl.CO.e_E(2)/k/T0+sum(nCOA_b)*Prcl.CO.e_E(3)/k/T0);
M(num_T,num_T)=2.5*(nm_b+sum(nCOa_b)+sum(nCOA_b)+nC2_b)...
                                                +1.5*(nac_b+nao_b+naAr_b);
AA = sparse(M);
    out=AA^(-1)*R;
end