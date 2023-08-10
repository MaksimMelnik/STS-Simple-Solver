function out = Discharge_DC_Hubner_air
% The main function for the macroparameters calculation in the afterglow
% of the electric discharge. Air plasma DC discharge in the tube.
% 0D simplification.
% Hubner's experiment conditions [1] accounting data and the system of
% equations from Pintassilgo et al [2].
% [1] M Hubner et al Meas. Sci. Technol. 23 (2012) 115602.
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.
% 5.06.2023 Maksim Melnik

%  todo:
% ions:
%   O2+
%   some O2+
%   O2+ reactions
% check URM
% check the paper
% e-
% turn on discharge
% add N(4S) +O+N2 → NO(X) + N2
% add N(4S) +O+O2 → NO(X) + O2
% ions:
%   O+
%   NO+
%   O-
% add back N2B
% add N2(B3Пg) formation reaction
% adding N2(B3Пg)
%   (R12) N2B + N2 -> N2A + N2
%   plot N2B + N2 -> N2A + N2
% send application
% fix n_N2B and Q_R7 N2B + O2 -> N2X + O + O
% - add all particles:
%   N2(B3Пg, B'3Σ−u, C3Пu, a'1Σ−u, a1Пg, w1Δu)
%   O2(a1Δg, b1Σ+g)
%   N(2D, 2P)
%   NO(A2Σ+, B2П)
%   NO2(X, A)
%   O3
%   N2+, N4+, O2+, O+, NO+, O−
%   e-
% - add all 15 reactions:
%   (R1)  e      + N2  → e+N∗2 → e+N(4S) + N(2D)
%   (R2)  e      + O2  → e+O2(A, C, c) → e+O(3P) + O(3P)
%   (R3)  e      + O2  → e+O2(B) → e+O(3P) + O(1D)
%   (R7)  N2(B)  + O2  → N2(X) + O + O
%   (R8)  N2(a') + O2  → N2(X) + O + O
%   (R9)  N2(a)  + O2  → N2(X) + O + O
%   (R10) N2(w)  + O2  → N2(X) + O + O
%   (R11) N2(A)  + O   → NO(X) + N(2D)
%   (R12) N2(B)  + N2  → N2(A) + N2
%   (R13) e      + N+2 → N(4S) + N(4S)
%   (R14) e+O+2 → O(3P) + O(3P)
%   (R15) e + NO+ → N(4S) + O(3P)
% - add 11 processes and corresponding Qin:
%   (1) Elastic collisions of electrons with N2 and O2
%   (2) Nitrogen and oxygen dissociation by electron
%   (4) VV N2-O2
%   (6) V–T energy exchanges in N2–N collisions involving multiquantum
%   (9) Diffusion of molecular and atomic metastable states to the wall
%   (10) Chemical reactions, Qchem
%   (11) Electron–ion recombination involving nitrogen or oxygen ions,Qe−i
% add vibrational excitation of O2
% add VT rates from V. Guerra works and fix the VT fluxes
% - N2-N    % first five transitions, same as i->i-1
% Consider the second Zeldovich reaction. It's not included in 
%   prof. Guerra's works, but it affects

warning("The present test case is unfinished")
warning(['Energy fluxes Q are not finnished, now only VT, VV, Diss, ' ...
    'Zeldovich, wall VT, wall rec, wall ET are included'])
warning('Thermal average velocity in R_VT_wall should be recheckerd')
warning('gamma_v in R_VT_wall is the same for each particle')
warning('Check energies in the wall recombination function')
warning('Zeldovich reactions are without electronic excitation')
warning("VV exchanges are with a crutch.")
warning('Find correct N2+ EM value.')
warning('Fix indexes in O2 + N -> NO + O, Kossyi1992')
disp('Started.')

tic                             % measuring computing time
    % constants
k = 1.380649e-23;               % Boltzmann constant, J/K
addpath('../src/')
load('../data/particles.mat', 'N2', 'O2', 'N', 'O', 'NO', 'N2p', 'O2p')
    % electronic excitation
N2.num_elex_levels = 1;         % N2(X1Σg+)
N2.num_elex_levels = 2;         % N2(X1Σg+, A3Σu+)
N2.num_vibr_levels(2) = 1;  N2.ev_0(2) = 0;  N2.ev_i{2} = 0;
% N2.num_elex_levels = 3;         % N2(X1Σg+, A3Σu+, B3Пg)
% N2.num_vibr_levels(3) = 1;  N2.ev_0(3) = 0;  N2.ev_i{3} = 0;
    % no electronic excitation
O2.num_elex_levels  = 1;         
O.num_elex_levels   = 1;
N.num_elex_levels   = 1;
NO.num_elex_levels  = 1;
N2p.num_elex_levels = 1;
O2p.num_elex_levels = 1;
    % no vibrational excitation
NO.num_vibr_levels(1) = 1;  NO.ev_0(1) = 0;  NO.ev_i{1} = 0;
O2.num_vibr_levels(1) = 1;  O2.ev_0(1) = 0;  O2.ev_i{1} = 0;

    % initial conditions
    % f_M_i are fractions of particle M at the moment i (i=0 is initial,
    %   i=3 is when the discharge is off. Fractions_3 are approximate.
init_c = [% p0, Pa; f_O2_0; f_NO_0; T0, K; T3, K; f_O_3; f_NO_3; f_N_3;
            133     0.2     0.008   300    440    1.2e-1 3.8e-3  2e-3 ...
        ... f_N2A_3 f_N2B_3 ion degree
            5.8e-5  7e-6    1e-6
        ];
for i_ini = 1           % choosing desired initial coonditions
 for i_U=3 % [2 3 4]    % choosing desired U dissociation parameter model
                        %   2 is for D/6k; 3 is for 3T; 4 is for inf
  for i_vibr = 3 % [1 2 3] % choosing vibrational energy exchange model
                        %   1 is for SSH; 2 is for FHO;
                        %   3 is for kinetics from V. Guerra works
   T0         = init_c(i_ini, 4);         % K
   n0         = init_c(i_ini, 1)/k/T0;    % m-3
   f_O2_0     = init_c(i_ini, 2);
   f_NO_0     = init_c(i_ini, 3);
   T3         = init_c(i_ini, 5) /T0;
   f_O_3      = init_c(i_ini, 6);
   f_N_3      = init_c(i_ini, 8);
   f_NO_3     = init_c(i_ini, 7);
   f_N2A_3    = init_c(i_ini, 9);
   f_N2B_3    = 0;
%    f_N2B_3 = init_c(i_ini, 10);
   ion_degree = init_c(i_ini, 11);
   n3      = init_c(i_ini, 1)/k/T3/T0; % m-3
   
   sigma0 = pi*N2.diameter^2;
   Delta = 1 / sqrt(2) / n0 / sigma0; % characteristic length, m
   t0    = 1 / (4 * n0 * N2.diameter^2 * sqrt(pi * k * T0 / N2.mass));

   num=0;
   Ps = {num, N2, O2, NO, N, O, N2p, O2p};
   index = cell(1, length(Ps));
   index{1}=0;
   for ind=2:length(Ps)
    num_v_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_v_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_v_states-1;
   end
   kinetics.Ps = Ps(2:end);
   kinetics.num_Ps = length(kinetics.Ps);
   
   Diss.Arrhenius='Park';
   Diss.rec=true;
   Diss.NEmodel='MT';
   switch i_U
    case 2
	 Diss.U='D/6k';
	case 3
	 Diss.U='3T';
	case 4
	 Diss.U='inf';
   end
   switch i_vibr
    case 1
	 model_VT = 'SSH';
	case 2
	 model_VT = 'FHO';
	case 3
	 model_VT = 'Guerra';
   end
   load('../data/reactions.mat', 'Reactions');
   ReactZel_1   = Reactions("N2 + O -> NO + N");
   React_N2A_O2 = Reactions("N2(A) + O2 -> N2(X) + O + O");
%    React_N2B_O2 = Reactions("N2(B) + O2 -> N2(X) + O + O");
   ReactZel_2 = Reactions("O2 + N -> NO + O");
%    Exch = [ReactZel_1("Kunova"), ReactZel_2("Kunova")];
%    Exch = [ReactZel_1("Kunova, NO(1)"), ReactZel_2("Kunova, NO(1)")];
        % V Guerra Zeldovich model
%    Exch = [ReactZel_1("Guerra95"), ReactZel_1("Guerra95_reverse"), ...
%                                                     ReactZel_2("Kunova")];
%    Exch = [ReactZel_1("Guerra95"), ReactZel_1("Guerra95_reverse")];
   Exch = [ReactZel_1("Guerra95"), ReactZel_1("Guerra95_reverse"), ...
       React_N2A_O2("Pintassilgo2009"), ReactZel_2("Kossyi1992")];
%    Exch = [ReactZel_1("Guerra95"), ReactZel_1("Guerra95_reverse"), ...
%        React_N2A_O2("Pintassilgo2009"), React_N2B_O2("Kossyi1992")];
   N2A_diff = Reactions("N2(A) + wall -> N2(X) + wall");
   ET_diff_c    = cell(1, kinetics.num_Ps);
                    % N2(X),          N2(A)
   ET_diff_c{1} = [Reactions("zero"), N2A_diff("Levron1978"), ...
       ...                                          N2(B)
                                                    Reactions("zero")];
   Reacs_keys = {'VT',     'VV',     'Exch', 'Wall', 'ET',    'Rec_wall'};
   reacs_val  = {model_VT, model_VT, Exch,   1,      ET_diff_c, 1};
   kinetics.num_eq = num;
   kinetics.reactions = containers.Map(Reacs_keys, reacs_val);
   kinetics.index = index(2:end);
   kinetics.n0 = n0;
   kinetics.T0 = T0;
   kinetics.Delta = Delta;
   kinetics.t0 = t0;
   kinetics.tube_R = 0.01;      % m
   kinetics.Tw     = 300;       % wall temperature, K
        % determine index numbers of molecules
   names=repmat("", length(kinetics.Ps), 1);
   serial_index=zeros(length(kinetics.Ps), 1);
   for i=1:length(kinetics.Ps)
        names(i)=string(kinetics.Ps{i}.name);
        serial_index(i)=i;
        IndexOfMolecules=containers.Map(names,serial_index);
   end
   kinetics.IndexOfMolecules=IndexOfMolecules;
%    xspan = [0.005 0.015]/t0;    % the Hubner experiment measurments time
   xspan = [0.005 0.2]/t0;      % from Pintassilgo 2014
%    xspan = [0.005 0.0052]/t0; % tests
   load('../data/for comparison/Hubner2012_and_Pintassilgo2014.mat' ...
                                                            ) %#ok<LOAD>
   i_vec = 0:30;
   N2_VDF = interp1(Pintassilgo2014_N2_VDF_post_DC(:, 1), ...
                        Pintassilgo2014_N2_VDF_post_DC(:, 2), i_vec, ...
                                        'spline', 'extrap');%#ok<USENS>
   N2_VDF(N2_VDF<0) = 0;
   n_N2 = N2.ev_i{1}*0;
   n_N2(i_vec+1) = N2_VDF/sum(N2_VDF);
   n_N2 = n_N2';
   Tv1 = N2.ev_i{1}(2)./(k*log(n_N2(1)./n_N2(2)));
   f_O2_3 = ((2-f_O_3-f_N_3)*(f_O2_0+f_NO_0/2) - f_O_3 - f_NO_3)/2;
   f_N2_3 = 1 - f_O2_3 - f_NO_3 - f_O_3 - f_N_3 - f_N2A_3 - f_N2B_3;
   n_N2 = n_N2 * f_N2_3 * (1 - ion_degree);
   n_O2 = density_f_exc(Tv1, f_O2_3, O2);
   n_NO = density_f_exc(Tv1, f_NO_3, NO);
       % N2(X,v), N2(A3Σu+), N2(B3Пg), O2(X), NO(X), N(X),  O(X),  
   y0 = [n_N2;    f_N2A_3;             n_O2;  n_NO;  f_N_3; f_O_3; ...
     ... N2+,               O2+
         f_N2_3*ion_degree; 0];
if N2.num_elex_levels == 3
   y0 = [n_N2;    f_N2A_3;   f_N2B_3;  n_O2;  n_NO;  f_N_3; f_O_3; ...
     ... N2+
         f_N2_3*ion_degree];
end
       % t3 correction, T
   y0 = [y0 * n3/n0;    T3];
%    y0 = [y0;            T3];
%    options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-20, ...
%                                     'NonNegative', 1:kinetics.num_eq+1); 
%    options_s = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, ...
%                                     'NonNegative', 1:kinetics.num_eq+1); 
%    options_s = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, ...
%                                     'NonNegative', 1:kinetics.num_eq+1); 
   options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+1); 
   [X, Y] = ode15s(@(t, y) ...
    Rpart_ODE_tube_DC_discharge_0D(t, y, kinetics), xspan, y0, options_s);

   t = X*t0;
   Y(:, 1:end-1) = Y(:, 1:end-1)*n0;
   Y(:, end)=Y(:, end)*T0;
   T=Y(:, end);
   Tv = N2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
   out.res=[t, Y, Tv];
   out.Y = Y;
   out.t = t;
   out.T = T;
   out.Tv = Tv;
   out.kinetics = kinetics;
  end
 end
end

%     figure  % T vs exp plot
% plot(Hubner_2012_T(:, 1), Hubner_2012_T(:, 2), 'sq', t*1e3, T, ...
%                         t*1e3, Tv, '-.', 'linewidth', 1.5) %#ok<USENS>
% % errorbar T Hubner +- 40 K
% legend('T_{exp}, Hubner 2012', 'T', 'Tv', 'location', 'best')
% xlabel('t, ms')
% xlim([-2 14])
% ylim([250 620])

rmpath('../src/')

addpath('../plots/')
Hubner12_Pintassilgo14_Air_DC_plots(out)
rmpath('../plots/')

toc
end