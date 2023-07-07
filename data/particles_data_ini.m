h = 6.626070041e-34;      % Plank constant, J*sec
N_a=6.02214076e23;        % Avogadro constant
c = 299792458;            % speed of light
k = 1.380649e-23;         % Boltzmann constant, J/K
warning('Only anharmonic oscillator parameters are used.')
addpath('../src/')


    % initialization of particles variables
C.name='C';
C.mass=1.994473440944e-26;          % kg
C.diameter=3.298e-10;   C.s_e=1;
C.m_mass=12.011;                    % molar mass
C.form_e=4.97472e-19;               % online table for kappa
C.form_e=1.180914618115280e-18;     % by Alena Kosareva
C.num_elex_levels=1;                % number of electronic levels
C.num_vibr_levels=1;                % actually atom has no vibr levels
C.e_E=0;                            % electronic excitation energy
C.fr_deg_c=3;                       % freedom degree for room temperature
C.EM=71.4;                          % Parameter ε/k (Lennard-Jones), К


O.name='O';
O.mass=2.6567628316576e-26;         % kg
O.m_mass=15.999;
O.diameter=2.75e-10;
O.s_e=5;                            % from DB, Capitelli
O.s_e=9;                            % staticsical weigth
O.m_mass=15.999;                    % molar mass
O.form_e= 4.098045681049634e-19;    % online table for kappa
O.form_e= 4.098111014876693e-19;    % by Alena Kosareva
O.num_elex_levels=1;                % number of electronic levels
O.num_vibr_levels=1;                % actually atom has no vibr levels
O.e_E=0;                            % electronic excitation energy
O.fr_deg_c=3;                       % freedom degree for room temperature
O.EM=80;                            % Parameter ε/k (Lennard-Jones), К
O.BMbeta=4.14;              % beta parameter of Born-Mayer potential, A^-1

N.name='N';
N.mass=2.32587E-26;
N.m_mass=14.0067;
N.s_e=4;
N.diameter=3.29800E-10;
N.num_elex_levels=1;                % number of electronic levels
N.num_vibr_levels=1;                % actually atom has no vibr levels
N.BMbeta=2.68;              % beta parameter of Born-Mayer potential, A^-1
N.EM=71.4;
N.fr_deg_c=3;
N.form_e=7.81808E-19;
N.e_E=0;                            % electronic excitation energy


CO.name='CO';
CO.mass=4.651236272601599e-26;            % molecular mass, kg
CO.m_mass=28.0104;                        % molar mass
CO.red_osc_mass=C.mass*O.mass/(C.mass+O.mass);    
% CO.diss_e=1.77716869535000e-18;         % dissociation energy
CO.diss_e=[89460 41121 24848]*100*h*c;    % dissociation energy for each 
                                          %           el. state, Krupenie
CO.diss_parts=["C", "O"];                 % dissociation products
% keys_Park={'CO', 'C', 'O', 'Ar', 'C2'}; %writing down Arrhenius law params
% val_Park_CO.A=2.3e20/N_a*1e-6;  val_Park_CO.N=-1;
% val_Park_C.A=3.4e20/N_a*1e-6;   val_Park_C.N=-1;
% val_Park_O.A=3.4e20/N_a*1e-6;   val_Park_O.N=-1;
% val_Park_Ar.A=2.3e19/N_a*1e-6;  val_Park_Ar.N=-1;
% val_Park_C2.A=2.3e20/N_a*1e-6;  val_Park_C2.N=-1;
% Map_Park=containers.Map(keys_Park, );
% CO.diss_Arrhenius=containers.Map(...
%     {'Park', 'Ibraguimova', 'McKenzie', 'Fairbairn'}, );
CO.form_e=-1.8349e-19;                    % DB
CO.form_e=-1.864225959573993e-19;         % by Alena Kosareva
CO.form_e_atoms_sum=C.form_e+O.form_e;
CO.num_elex_levels=3;                     % number of electronical levels
CO.num_vibr_levels=[68 34 18];            % number of vibrational levels
CO.diameter=3.65e-10;                     % molecular diameter
    % фактор симметрии (1 для гетероядерных, 2 для гомоядерных)
CO.sigma=1;                         
CO.Be=[193.128087 169.124 161.15];        % спектроскопическая постоянная
CO.moment_I=1.4560E-46;
CO.mltpl_atoms_mass=C.mass*O.mass;        % произведение масс атомов
CO.mltpl_atoms_s_e=C.s_e*O.s_e;           % произведение статвесов атомов
CO.s_e=[1 6 2];                           % electronic statistical weight2
    % спектроскопические постоянные
CO.we = [216981.358 174341 151824];       % DB, м-1
CO.wexe=[1328.831 1436 1940]; 
CO.weye=[1.0511 -4.5 76.6];         
CO.e_E=[0  4868740*h*c  6507460*h*c];     % electronic excitation 
                                          %             energy, Krupenie
    % vibrational levels energy
e_i=levels_e_ex(CO, 1);
CO.ev_0(1)=e_i(1);  CO.ev_i{1}=e_i-e_i(1);  % the new e_vibr standart
e_i=levels_e_ex(CO, 2);
CO.ev_0(2)=e_i(1);  CO.ev_i{2}=e_i-e_i(1);  % the new e_vibr standart
e_i=levels_e_ex(CO, 3);
CO.ev_0(3)=e_i(1);  CO.ev_i{3}=e_i-e_i(1);  % the new e_vibr standart
CO.C_VE=[0 13.58e-2 5.49e-2];                   % VE-exchange parameters
CO.S_VE=[0 5e-13 9.94e-13]/1e6;                 % in SI
CO.mA_mAB=0.42880501528003884;  CO.mB_mAB=0.5711949847199612;  % FHO param
CO.fr_deg_c=5;                      % freedom degree for room temperature
CO.EM=98.1;                         % Parameter ε/k (Lennard-Jones), К
CO.r_e=1.128323e-10;                % internuclear distance, r_e, m


C2.name='C2';
C2.mass=3.988946881888e-26;         C2.diameter=3.621e-10; 
C2.m_mass=2*C.m_mass;                  % molar mass
C2.diss_e=9.949436628352846e-19;
C2.diss_parts=["C", "C"];
C2.ev_i{1}=0;      C2.ev_0=0;
C2.e_E=0;
C2.Be=181.98399999999998;   C2.sigma=2;
C2.s_e=1;   C2.e_E=0;   
C2.num_elex_levels=1;   
C2.num_vibr_levels=1;
C2.fr_deg_c=5;                         % freedom degree for room temperature
C2.mltpl_atoms_mass=C.mass*C.mass;     % произведение масс атомов
C2.mltpl_atoms_s_e=C.s_e*C.s_e;        % произведение статвесов атомов
C2.form_e=0;                           % Formation energy
C2.form_e_atoms_sum=C.form_e+C.form_e;
C2.CO_C_C2_O_e=58000*k;
C2.EM=97.53;                           % DB


Ar.name='Ar';
Ar.mass=6.633521356992e-26; Ar.diameter=3.33e-10;
Ar.m_mass=39.948;                   % molar mass
Ar.num_elex_levels=1;               % number of electronic levels
Ar.num_vibr_levels=1;               % actually atom has no vibr levels
Ar.e_E=0;                           % electronic excitation energy
Ar.fr_deg_c=3;                      % freedom degree for room temperature
Ar.form_e=0;                        % Formation energy
Ar.EM=136.5;                        % DB


O2.name='O2';                       % data from DB work-v5
O2.mass=5.31353E-26;                % kg
O2.m_mass=31.9988;                  % molar mass
O2.red_osc_mass=0.5*O.mass;
O2.diameter=3.5155E-10;             % m
O2.sigma=2;
O2.Be=[     143.768,    142.64,     140.03699999999998, 91.55, 96.0, ...
            91.06,      53.0,       81.89999999999999, 49.0, 71.17, ...
            57.06,      81.89999999999999, 78.8, 54.7,      52.59, ...
            173.0,      170.3,      146.38,     81.10000000000001]; % m-1
O2.form_e=0;
O2.form_e_atoms_sum=2*O.form_e;
O2.num_elex_levels=1;               % number of electronical levels
O2.num_vibr_levels=36;              % number of vibrational levels
O2.diss_e=[ 8.19609E-19 6.63036E-19 5.73147E-19 1.78601E-19 1.46163E-19...
            1.32119E-19 5.50047E-20 1.59798E-19 7.56180E-20 5.23826E-20...
            2.35592E-20 1.78005E-19 1.44991E-19 7.00222E-20 3.14057E-20...
            4.39024E-19 4.22517E-19 3.47191E-19 2.64853E-19];   % J

%Data from Park 1994 
keys={'Ar', 'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2'};
val_A=[2e21, 1e22, 1e22, 1e22, 2e21, 2e21, 2e21, 2e21, 2e21, ...
                                                    2e21, 2e21]/N_a*1e-6;
val_n = num2cell(val_A*0 - 1.5);
val_A = num2cell(val_A);
O2.diss_Arrhenius_A=containers.Map(keys, val_A);
O2.diss_Arrhenius_n=containers.Map(keys, val_n);
O2.diss_parts=["O", "O"];
O2.mltpl_atoms_mass=O.mass^2;
O2.mltpl_atoms_s_e=O.s_e(1)^2;
O2.s_e=[    3,          2,          1,          1,          6, ...
            3,          10,         3,          3,          6, ...
            6,          8,          8,          2,          2, ...
            3,          1,          3,          1];
O2.we=[     158019.00   148350.00   143277.00   79429.00    85000.00 ...
            79907.00    20000.00    70931.00    53700.00    20000.00 ...
            20000.00    162640.00   70560.00    49950.00    20000.00 ...
            195700.00   192700.00   254700.00   79240.00];      % m-1
O2.wexe=[   1198.00     1290.00     1400.00     1273.60     2000.00 ...
            1216.00     0.00        1065.00     1373.00     0.00 ...
            0.00        16370.00    960.00      1350.00     0.00 ...
            1970.00     1900.00     0.00        770.00];        % m-1
O2.weye=[   4.747       0.000       0.000       -24.440     0.000 ...
            -55.000     0.000       -13.900     0.000       0.000 ...
            0.000       0.000       0.000       0.000       0.000 ...
            0.000       0.000       0.000       0.000];         % m-1
O2.e_E=0;
O2.r_e=1.20752e-10;                 % internuclear distance, m
O2=ev_i_ini(O2);                    % vibr energy
O2.fr_deg_c=5;                      % freedom degree for room temperature
O2.EM=107.4;                        % Parameter ε/k (Lennard-Jones), К
O2.BMbeta=3.964;            % beta parameter of Born-Mayer potential, A^-1


N2.name='N2';                       % data from DB work-v5
N2.mass=4.65173E-26;                % kg
N2.m_mass=28.0134;
N2.red_osc_mass=0.5*N.mass;
N2.diameter=3.4039E-10;             % m
N2.sigma=2;
N2.Be=[199.824      145.460     163.745     147.000...
147.330     147.990     161.690     149.800...
92.100      92.800      182.473];
N2.form_e=0;
N2.form_e_atoms_sum=2*N.form_e;
N2.num_elex_levels=1;               % number of electronical levels
N2.num_vibr_levels=47;              % number of vibrational levels
N2.diss_e=[ 1.56362E-18 5.89862E-19 7.84447E-19 7.80693E-19 8.43286E-19...
            9.96818E-19 9.74372E-19 9.18473E-19 6.84927E-20 2.22661E-19...
            1.98088E-19];           % J
%Data from Park 1994 
keys={'Ar', 'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2'};
val_A=[7e21, 3e22, 3e22, 3e22, 7e21, 7e21, 7e21, 7e21, 7e21, ...
                                                    7e21, 7e21]/N_a*1e-6;
val_n = num2cell(val_A*0 - 1.6);
val_A = num2cell(val_A);
N2.diss_Arrhenius_A=containers.Map(keys, val_A);
N2.diss_Arrhenius_n=containers.Map(keys, val_n);
N2.diss_parts=["N", "N"];
N2.mltpl_atoms_mass=N.mass^2;
N2.mltpl_atoms_s_e=N.s_e(1)^2;
N2.s_e=[1,      3,       6,      6,...
3,      1,      2,      2,...
5,      6,      6];
N2.we=[     235857      146064      173339      150140      151688 ...
            153025      169420      155926      66700       74249 ...
            204718];                % m-1
N2.wexe=[   1432.40     1387.00     1412.20     1160.00     1218.00 ...
            1207.00     1394.90     1163.00     0.00        1185.00 ...
            2844.50];               % m-1
N2.weye=[   -0.226      1.030       -5.690      0.000       4.186 ...
            4.129       0.794       0.000       0.000       0.000 ...
            208.833];               % m-1
N2.e_E=0;
N2.r_e=1.09768E-10;                 % internuclear distance, m
N2=ev_i_ini(N2);                    % vibr energy
N2.fr_deg_c=5;                      % freedom degree for room temperature
N2.EM=97.53;                        % Parameter ε/k (Lennard-Jones), К
N2.BMbeta=2.573;            % beta parameter of Born-Mayer potential, A^-1



NO.name='NO';
NO.mass=4.98263E-26;
NO.m_mass=30.0061;                  % molar mass
NO.red_osc_mass=O.mass*N.mass/(O.mass+N.mass);
NO.diameter=3.4061E-10;
NO.sigma=1;
NO.Be=[167.195  112.750  199.650    112.440...
133.580     200.000     200.260     112.750...
133.200     198.630     198.200     200.300...
125.230     113.200     203.400]; % m-1
NO.form_e=1.50711E-19;
NO.form_e_atoms_sum=N.form_e+O.form_e;
NO.num_elex_levels=1;               % number of electronic levels
NO.num_vibr_levels=38;              
NO.diss_e=[1.04090E-18      2.88770E-19     5.68183E-19     5.29124E-19...
4.20809E-19     4.04897E-19     3.87059E-19     3.76074E-19...
2.42446E-19     2.37182E-19     2.13920E-19     2.00552E-19...
1.91791E-19     3.80285E-19     3.59666E-19];   % J
%Data from Park 1994 
keys={'Ar', 'C', 'N', 'O', 'C2', 'N2', 'O2', 'CN', 'CO', 'NO', 'CO2'};
val_A=[5e15, 1.1e17, 1.1e17, 1.1e17, 5e15, 5e15, 5e15, 5e15, 5e15, ...
                                                    1.1e17, 1.1e17]/N_a*1e-6;
val_n = num2cell(val_A*0 - 0);
val_A = num2cell(val_A);
NO.diss_Arrhenius_A=containers.Map(keys, val_A);
NO.diss_Arrhenius_n=containers.Map(keys, val_n);
NO.diss_parts=["N", "O"];
NO.mltpl_atoms_mass=O.mass*N.mass;
NO.mltpl_atoms_s_e=O.s_e*N.s_e;
NO.s_e=[4 ,      8 ,      2 ,      4,...
4,       4,       2 ,      4 ,      4,...
2,       4,       2 ,      2 ,      4 ,      4];
NO.we=[190420.00        101900.00       237471.00       103720.00...
120600.00       239500.00       232390.00       100440.00...
121740.00       237530.00       239400.00       233940.00...
108554.00       95200.00        243830.00];      % m-1
NO.wexe=[1407.50        1280.00     1610.60     777.00...
1500.00     1500.00     2288.50     1100.00...
1561.00     1640.00     2000.00     0.00...
1108.30     1128.00     4838.00];        % m-1
NO.weye=[1.100      0.000       -4.650      10.000...
0.000       0.000       75.000      0.000...
0.000       0.000       0.000       0.000...
-14.390     0.000       0.000];         % m-1
NO.e_E=0;
NO.r_e=1.15*1e-10;                 % internuclear distance, m
NO=ev_i_ini(NO);                    % vibr energy
NO.fr_deg_c=5;                      % freedom degree for room temperature
NO.EM=119;                      % Parameter ε/k (Lennard-Jones), К
NO.BMbeta=3.303;            % beta parameter of Born-Mayer potential, A^-1




    % объявляем столкновения
    Coll_CO_CO.coll_diameter=CO.diameter;
% Park
    Coll_CO_CO.ArrA(1)=2.3e20/N_a*1e-6;     Coll_CO_CO.ArrN(1)=-1;
    Coll_CO_C.ArrA(1) =3.4e20/N_a*1e-6;     Coll_CO_C.ArrN(1) =-1;
    Coll_CO_O.ArrA(1) =3.4e20/N_a*1e-6;     Coll_CO_O.ArrN(1) =-1;
    Coll_CO_Ar.ArrA(1)=2.3e19/N_a*1e-6;     Coll_CO_Ar.ArrN(1)=-1;
    Coll_CO_C2.ArrA(1)=2.3e20/N_a*1e-6;     Coll_CO_C2.ArrN(1)=-1;
    Coll_C2.ArrA(1)   =3.7e14/N_a*1e-6;     Coll_C2.ArrN(1)   = 0;
    Coll_CO_C__C2_O.ArrA(1)=2e17/N_a*1e-6;  Coll_CO_C__C2_O.ArrN(1)=-1;
    Coll_NO_O__N_O2.ArrA(1)=8.4e12/N_a*1e-6; Coll_NO_O__N_O2.ArrN(1)=0;
    Coll_N2_O__NO_N.ArrA(1)=6.4e17/N_a*1e-6; Coll_N2_O__NO_N.ArrN(1)=-1;
% Ibraguimova
    Coll_CO_CO.ArrA(2)=2.8e21/N_a*1e-6;     Coll_CO_CO.ArrN(2)=-1.39;
    Coll_CO_C.ArrA(2) =1.4e22/N_a*1e-6;     Coll_CO_C.ArrN(2) =-1.39;
    Coll_CO_O.ArrA(2) =2.1e22/N_a*1e-6;     Coll_CO_O.ArrN(2) =-1.39;
    Coll_CO_Ar.ArrA(2)=1.4e21/N_a*1e-6;     Coll_CO_Ar.ArrN(2)=-1.39;
% МакКензи
    Coll_CO_CO.ArrA(3)=0;	Coll_CO_CO.ArrN(3)=0;
% Fairbairn diss coefficients
    CO_CO_krec=2e-34/1e12;  % cm6 s-1, Fairbairn
    Tfb=7600;               % Temperature in Fairbairn experiment
    Zvibr=sum(exp(-CO.ev_i{1}'/(k*Tfb)));
    Io2=1.95790474021232e-46;
    Zrot=8*pi^2*Io2*k*Tfb/(2*h^2);
    Zint=Zrot.*Zvibr;
    kdr_fb=(CO.mass/CO.mltpl_atoms_mass)^1.5*h^3*...
                (2*pi*k*Tfb).^(-1.5).*Zint.*exp(CO.diss_e(1)./(k*Tfb));
    Coll_CO_CO.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_CO.ArrN(4)=0;
    Coll_CO_C.ArrA(4) =CO_CO_krec/kdr_fb;    Coll_CO_C.ArrN(4) =0;
    Coll_CO_O.ArrA(4) =CO_CO_krec/kdr_fb;    Coll_CO_O.ArrN(4) =0;
    Coll_CO_Ar.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_Ar.ArrN(4)=0;
    Coll_CO_C2.ArrA(4)=CO_CO_krec/kdr_fb;    Coll_CO_C2.ArrN(4)=0;
    Coll_C2.ArrA(4)   =6e-10/1e6;            Coll_C2.ArrN(4)   =0;
    Coll_CO_C__C2_O.ArrA(4)=6e-10/1e6;       Coll_CO_C__C2_O.ArrN(4)=0;
% Appleton diss coefficients with global diss energy
    Coll_CO_CO.ArrA(5)=2.8e-3/1e6;      Coll_CO_CO.ArrN(5)=-2.86;
    Coll_CO_C.ArrA(5) =2.8e-3/1e6;      Coll_CO_C.ArrN(5) =-2.86;
    Coll_CO_O.ArrA(5) =2.8e-3/1e6;      Coll_CO_O.ArrN(5) =-2.86;
    Coll_CO_Ar.ArrA(5)=2.8e-3/1e6;      Coll_CO_Ar.ArrN(5)=-2.86;
    Coll_CO_C2.ArrA(5)=2.8e-3/1e6;      Coll_CO_C2.ArrN(5)=-2.86;
    Coll_C2.ArrA(5)=3.7e14/N_a*1e-6;    Coll_C2.ArrN(5)   =0;
% Appleton diss coefficients with own diss energy
% CO.diss_e(1)=98600*k;
%     Coll_CO_CO.ArrA(6)=4.4e-10/1e6;     Coll_CO_CO.ArrN(6)=0;
%     Coll_CO_C.ArrA(6) =4.4e-10/1e6;     Coll_CO_C.ArrN(6) =0;
%     Coll_CO_O.ArrA(6) =6.86e-9/1e6;     Coll_CO_O.ArrN(6) =0;
%     Coll_CO_Ar.ArrA(6)=4.4e-10/1e6;     Coll_CO_Ar.ArrN(6)=0;
%     Coll_CO_C2.ArrA(6)=4.4e-10/1e6;     Coll_CO_C2.ArrN(6)=0;
%     Coll_C2.ArrA(6)=3.7e14/N_a*1e-6;    Coll_C2.ArrN(6)=   0;
    
% Mick diss coefficients=32 
    Coll_CO_CO.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_CO.ArrN(7)=-3.1;
    Coll_CO_C.ArrA(7) =4.3e27/N_a*1e-6;     Coll_CO_C.ArrN(7) =-3.1;
    Coll_CO_O.ArrA(7) =4.3e27/N_a*1e-6;     Coll_CO_O.ArrN(7) =-3.1;
    Coll_CO_Ar.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_Ar.ArrN(7)=-3.1;
    Coll_CO_C2.ArrA(7)=4.3e27/N_a*1e-6;     Coll_CO_C2.ArrN(7)=-3.1;
    Coll_C2.ArrA(7)   =3.72e14/N_a*1e-6;    Coll_C2.ArrN(7)   =0;
    Coll_CO_C__C2_O.ArrA(7)=6e-10/1e6;       Coll_CO_C__C2_O.ArrN(7)=0; 

save particles.mat C O N CO C2 Ar O2 N2 NO
%save CO_C_O_Ar_C2.mat CO C O Ar C2
%save O2_O O2 O

addpath('../src/')

function out=ev_i_ini(M)
% automatic vibrational energy calculation
    for i_e=1:M.num_elex_levels
        e_i=levels_e_ex(M, i_e);
        M.ev_0(i_e)=e_i(1);
        M.ev_i{i_e}=e_i-e_i(1);
    end
    out=M;
end