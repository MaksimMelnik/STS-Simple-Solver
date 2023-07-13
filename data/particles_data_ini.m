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
O.form_e= 4.098045681049634e-19;    % DBparticles
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
    % staticsical weigth
N.s_e=[4 10 6 12 6 12 2 20 12 4 10 6 10 12 6 6 28 12 14 20 10 2 20 12 10 4 6 12 6 6 28 12 14 20 126 10 10 14 2 6 20 12 4 10 6 12 6 6 28 12 14 20 10 126 162 6 2 20 12 4 10 12 6 6 28 12 14 20 10 126 162 6 198 2 20 12 4 10 12 6 6 12 28 14 20 10 126 162 6 198 2 234 12 0 4 6 10 12 6 28 14 12 10 20 126 6 162 198 2 234 270 12 20 4 6 20 6 10 14 12 28 10 12 126 6 162 198 2 234 270 306 12 20 4 6 12 6 14 20 10 28 10 126 6 12 162 198 2 2 234 270 306 342 12 6 20 12 4 6 14 20 10 28 10 126 12 6 162 198 12 2 6 234 270 306 342 12 20 4 378 6 14 20 10 28 12 126 6 12 162 6 198 2 234 270 306 12 342 20 6 414 4 10 378 14 20 10 28 12 12 126 6 6 162 198 2 234 270 306 12 342 6 20 414 4 450 378 14 20 10 10 28 12 12 6 126 6 162 198 2 234 270 6 306 342 12 20 486 378 450 414 4 14 20 28 10 10 12 12 6 126 6 162 198 234 6 2 270 306 342 12 20 522 378 414 486 4 450 14 20 12 28 10 12 6 6 126 10 162 198 6 2 234 270 306 342 12 20 4 414 522 378 14 558 486 450 20 12 28 10 6 12 6 126 162 198 6 2 10 234 270 306 12 342 20 450 486 594 522 4 14 414 378 558 20 12 28 10 6 12 126 6 162 6 198 2 234 270 306 342 12 10 14 20 558 522 486 414 450 630 594 378 4 20 28 10 12 126 6 162 6 198 234 2 270 342 12 306 14 20 414 630 522 378 450 594 4 486 558 666 10 20 28 10 12 126 162 198 234 306 270 342 378 414 558 486 522 666 630 702 450 594 10 7938 10 8712 9522 10368 11250 12168 13122 14112 15138 16200 17298 18432 19602 20808 22050 23328 24642 25992 27378 28800 30258 31752 33282 34848 36450 38088 39762 41472 43218 45000 46818 48672 50562 52488 54450 56448 58482 60552 62658 64800 66978 69192 71442 73728 76050 78408 80802 83232 85698 88200 90738 95922 93312 98568 101250 103968 109512 106722 112338 118098 115200 124002 127008 121032 133128 130050 136242 139392 142578 145800 149058 152352 155682 159048 162450 165888 169362 172872 176418 180000 183618 187272 190962 194688 198450 202248 206082 209952 213858 217800 221778 225792 229842 238050 233928 259200 242208 246402 250632 254898 276768 281250 272322 263538 267912 294912 299538 285768 290322 323208 328050 332928 318402 304200 308898 313632 357858 352800 347778 337842 342792 399618 388962 394272 405000 421362 415872 410418 383688 438048 432450 426888 362952 378450 373248 368082 563922 557568 443682 455058 449352 532512 526338 538722 551250 544968 460800 502002 496008 490050 520200 514098 508032 478242 472392 466578 484128 656658 670482 663552 642978 636192 649800 677448 705672 720000 712818 684450 691488 698562 616050 583200 609408 629442 622728 570312 576738 602802 589698 596232 10 14 18 10 6 2 10 14 6 10 14 18 10 70 6 2 10 14 6 10 14 18 10 70 6 2 10 14 6 90 10 14 10 18 70 6 2 10 14 110 6 90 10 14 18 10 70 6 2 10 14 110 6 130 90 10 14 18 10 70 6 2 10 10 14 110 6 130 150 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 190 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 210 190 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 230 190 210 90 10 14 18 10 70 6 2 10 14 110 6 130 150 170 210 190 230 250 90 10 14 18 10 70 6 2 10 14 110 6 130 150 170 190 210 230 270 250 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 190 270 230 290 250 210 90 10 14 10 18 70 6 2 10 14 110 6 130 150 170 290 250 210 270 190 310 230 90 10 14 70 10 18 6 2 10 14 110 6 130 150 170 310 230 290 190 210 330 270 90 10 14 18 70 10 6 2 10 14 110 6 130 150 170 270 330 190 250 350 210 230 290 310 90 10 14 10 70 18 6 2 10 14 110 6 130 150 170 230 190 210 250 270 330 370 350 310 290 90 14 10 6 70 18 2 110 130 170 150 290 190 270 310 330 210 370 390 350 230 250 90 4410 4840 5290 5760 6250 6760 250 7290 7840 8410 9000 9610 10240 10890 11560 12250 12960 13690 14440 15210 16000 16810 17640 18490 19360 20250 21160 22090 23040 24010 25000 26010 27040 28090 29160 30250 31360 32490 33640 34810 36000 38440 37210 39690 40960 42250 43560 44890 46240 47610 50410 49000 51840 53290 54760 56250 57760 59290 60840 62410 64000 65610 67240 68890 70560 72250 73960 75690 77440 79210 81000 82810 84640 86490 88360 90250 92160 94090 96040 98010 100000 102010 104040 106090 108160 110250 112360 114490 116640 118810 121000 123210 125440 127690 129960 132250 134560 136890 139240 141610 144000 158760 156250 153760 151290 148840 146410 161290 163840 166410 169000 171610 174240 176890 179560 182250 184960 187690 190440 193210 196000 198810 243360 240250 237160 234090 231040 228010 225000 222010 219040 216090 213160 210250 207360 204490 201640 289000 299290 292410 313290 309760 256000 252810 246490 295840 249640 302760 278890 262440 282240 275560 272250 285610 306250 265690 259210 268960 327610 331240 388090 357210 324000 334890 349690 316840 320410 384160 368640 361000 372490 376360 364810 400000 392040 380250 396010 353440 342250 338560 345960 2 18 2 2 12 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 882 968 1058 1152 1250 1352 1458 1568 1682 1800 1922 2048 2178 2312 2450 2592 2738 2888 3042 3200 3362 3528 3698 3872 4050 4232 4418 4608 4802 5000 5202 5408 5618 5832 6050 6272 6498 6728 6962 7200 7442 7688 7938 8192 8450 8712 8978 9248 9522 9800 10082 10368 10658 10952 11250 11552 11858 12168 12482 12800 13122 13448 13778 14112 14450 14792 15138 15488 15842 16200 16562 16928 17298 17672 18050 18432 18818 19602 20000 20808 19208 20402 21218 21632 22050 22472 22898 23328 23762 24200 24642 25088 25538 25992 26450 26912 27378 27848 28322 28800 29282 29768 30258 30752 31250 31752 32258 32768 33282 33800 34322 34848 35378 35912 36450 36992 37538 38088 38642 39200 39762 40328 40898 41472 42050 42632 43218 43808 44402 45000 45602 46208 46818 47432 48050 48672 51200 50562 57122 55778 57800 58482 56448 61250 60552 59858 59168 62658 49298 51842 49928 61952 54450 52488 53138 55112 53792 79202 64082 77618 66978 78408 73728 76050 63368 75272 76832 74498 67712 71442 69938 72200 69192 72962 70688 65522 64800 68450 80000 66248 10 18 12 10 18 12 6 10 18 12 10 18 12 10 12 18 10 12 18 10 12 18 10 12 18 12 10 18 12 10 18 12 10 12 18 10 12 18 12 10 18 12 10 18 12 10 18 12 10 18 10 18 4410 4840 5290 5760 6250 6760 7290 7840 8410 9000 9610 10240 10890 11560 12250 12960 13690 14440 15210 16000 16810 17640 18490 19360 20250 21160 22090 23040 24010 25000 26010 27040 28090 29160 30250 31360 32490 33640 34810 36000 38440 37210 39690 40960 42250 43560 46240 44890 47610 50410 49000 53290 51840 54760 56250 57760 60840 62410 59290 65610 64000 68890 70560 67240 73960 72250 81000 75690 77440 79210 86490 82810 84640 94090 88360 90250 92160 102010 100000 104040 96040 98010 108160 106090 112360 110250 114490 123210 116640 118810 125440 121000 127690 141610 136890 139240 144000 132250 129960 134560 161290 166410 158760 156250 151290 153760 163840 148840 146410 171610 193210 169000 179560 174240 198810 184960 187690 182250 196000 190440 176890 201640 240250 243360 204490 237160 219040 207360 225000 210250 213160 216090 234090 222010 231040 228010 302760 299290 306250 275560 246490 256000 252810 313290 285610 292410 289000 278890 309760 295840 282240 272250 262440 265690 259210 249640 268960 345960 320410 353440 327610 349690 316840 324000 331240 342250 334890 338560 388090 357210 361000 364810 368640 372490 376360 400000 380250 384160 392040 396010 6];                            
N.form_e = 7.81808E-19;             % DBparticles
N.num_elex_levels=1;                % number of electronic levels
N.num_vibr_levels=1;                % actually atom has no vibr levels
N.e_E=0;                            % electronic excitation energy
N.fr_deg_c=3;                       % freedom degree for room temperature
N.EM=71.4;                          % Parameter ε/k (Lennard-Jones), К
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
O2.num_elex_levels = length(O2.Be); % number of electronical levels
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
    % electronic excitation energy, J
O2.e_E = [  0.00000E+00 1.57289E-19 2.62114E-19 6.56665E-19 6.89098E-19...
            7.03158E-19 7.80256E-19 9.89117E-19 1.07330E-18 1.10295E-18...
            1.13309E-18 1.34762E-18 1.37309E-18 1.44697E-18 1.48717E-18...
            1.49500E-18 1.51151E-18 1.58683E-18 1.63218E-18];
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
N2.Be=[     199.824     145.460     163.745     147.000     147.330 ...
            147.990     161.690     149.800     92.100      92.800  ...
            182.473];
N2.form_e=0;
N2.form_e_atoms_sum=2*N.form_e;
N2.num_elex_levels = length(N2.Be); % number of electronical levels
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
N2.s_e=[    1           3           6           6           3 ...
            1           2           2           5           6 ...
            6];
N2.we=[     235857      146064      173339      150140      151688 ...
            153025      169420      155926      66700       74249 ...
            204718];                % m-1
N2.wexe=[   1432.40     1387.00     1412.20     1160.00     1218.00 ...
            1207.00     1394.90     1163.00     0.00        1185.00 ...
            2844.50];               % m-1
N2.weye=[   -0.226      1.030       -5.690      0.000       4.186 ...
            4.129       0.794       0.000       0.000       0.000 ...
            208.833];               % m-1
N2.e_E = [  0.00000E+00 9.97267E-19 1.18431E-18 1.18805E-18 1.31647E-18...
            1.35382E-18 1.37627E-18 1.43217E-18 1.51836E-18 1.74609E-18...
            1.77064E-18];            % electronic excitation energy, J
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
NO.num_elex_levels = length(NO.Be);          % number of electronic levels
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
NO.mltpl_atoms_s_e = O.s_e(1) * N.s_e(1);
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
    % electronic excitation energy, J
NO.e_E = [  0.00000E+00 7.70880E-19 8.73355E-19 9.12414E-19 9.67002E-19...
            1.03545E-18 1.05450E-18 1.06547E-18 1.19910E-18 1.20436E-18...
            1.22762E-18 1.24100E-18 1.24973E-18 1.25226E-18 1.27286E-18];
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

rmpath('../src/')

function out=ev_i_ini(M)
% automatic vibrational energy calculation
    for i_e=1:M.num_elex_levels
        M.num_vibr_levels(i_e) = 100;
        e_i=levels_e_ex(M, i_e);
        temp = find( e_i(2:end) > e_i(1:end-1), 1, 'last') + 1;
        num_lvls = find(e_i(1:temp) < M.diss_e(i_e), 1, 'last');
        M.num_vibr_levels(i_e) = num_lvls;
        e_i = e_i(1:num_lvls);
        M.ev_0(i_e)=e_i(1);
        M.ev_i{i_e}=e_i-e_i(1);
    end
    out=M;
end