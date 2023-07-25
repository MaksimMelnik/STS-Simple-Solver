function out=Air_discharge_Hubner_0D
% The main function for the macroparameters calculation in the afterglow
% of the electric discharge.
% Hubner's experiment conditions [1] accounting data Pintassilgo et al [2]
% [1] M Hubner et al Meas. Sci. Technol. 23 (2012) 115602.
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.
% 15.02.2023 Maksim Melnik

tic % measuring computing time
    % constants
k = 1.380649e-23;             % Boltzmann constant, J/K
addpath('../src/')
load('../data/particles.mat', 'N2', 'O2', 'N', 'O', 'NO')
O2.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;
N2.num_elex_levels=1;
N.num_elex_levels=1;
NO.num_elex_levels=1;

    % initial conditions
    % f_M_i are fractions of particle M at the moment i (i=0 is initial,
    %   i=3 is when the discharge is off. Fractions_3 are approximate.
init_c=[% p0, Pa; f_O2_0; f_NO_0; T0, K; T3, K; f_O_3; f_NO_3; f_N_3 
          133     0.2     0.008   300    440    1.2e-1 3.8e-3  2e-3
        ];
for i_ini=1             % choosing desired initial coonditions
 for i_U=3 % [2 3 4]    % choosing desired U dissociation parameter model
                        %   2 is for D/6k; 3 is for 3T; 4 is for inf
  for i_vibr=1 % [1 2]  % choosing vibrational energy exchange model
                        %   1 is for SSH; 2 is for FHO
   T0     = init_c(i_ini, 4);         % K
   n0     = init_c(i_ini, 1)/k/T0;    % m-3
   f_O2_0 = init_c(i_ini, 2);
   f_NO_0 = init_c(i_ini, 3);
   T3     = init_c(i_ini, 5) /T0;
   f_O_3  = init_c(i_ini, 6);
   f_N_3  = init_c(i_ini, 8);
   f_NO_3 = init_c(i_ini, 7);
   
   sigma0 = pi*N2.diameter^2;
   Delta = 1 / sqrt(2) / n0 / sigma0; % characteristic length, m
   t0    = 1 / (4 * n0 * N2.diameter^2 * sqrt(pi * k * T0 / N2.mass));

   num=0;
   index{1}=0;
   Ps={num, N2, O2, NO, N, O};
   for ind=2:length(Ps)
    num_v_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_v_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_v_states-1;
   end
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
	 model_VT='SSH';
	case 2
	 model_VT='FHO';
   end
   Exch=1;
   Reacs_keys={'VT'};
   reacs_val={model_VT};
   Reacs_keys={'Diss', 'VT'};
   reacs_val={Diss, model_VT};
   Reacs_keys={'Diss', 'VT', 'VV'};
   reacs_val={Diss, model_VT, model_VT};
   Reacs_keys={'Diss', 'VT', 'VV', 'Exch'};
   reacs_val={Diss, model_VT, model_VT, Exch};
   kinetics.Ps=Ps(2:end);
   kinetics.num_Ps=length(kinetics.Ps);
   kinetics.num_eq=num;
   kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
   kinetics.index=index(2:end);
   kinetics.n0=n0;
   kinetics.T0=T0;
   kinetics.Delta=Delta;
   kinetics.t0=t0;
        %determine index numbers of molecules
   names=repmat("", length(kinetics.Ps), 1);
   serial_index=zeros(length(kinetics.Ps), 1);
   for i=1:length(kinetics.Ps)
        names(i)=string(kinetics.Ps{i}.name);
        serial_index(i)=i;
        IndexOfMolecules=containers.Map(names,serial_index);
   end
   kinetics.IndexOfMolecules=IndexOfMolecules;
   xspan=[0.005 0.015]/t0;      % the Hubner experiment measurments time
   load('../data/for comparison/Hubner2012_and_Pintassilgo2014.mat' ...
                                                            ) %#ok<LOAD>
   i_vec=0:30;
   N2_VDF = interp1(Pintassilgo2014_N2_VDF_post_DC(:, 1), ...
                                Pintassilgo2014_N2_VDF_post_DC(:, 2), ...
                                i_vec, 'spline', 'extrap'); %#ok<USENS>
   N2_VDF(N2_VDF<0) = 0;
   n_N2 = N2.ev_i{1}*0;
   n_N2(i_vec+1) = N2_VDF/sum(N2_VDF);
   n_N2 = n_N2';
   Tv1 = N2.ev_i{1}(2)./(k*log(n_N2(1)./n_N2(2)));
   f_O2_3 = ((2-f_O_3-f_N_3)*(f_O2_0+f_NO_0/2) - f_O_3 - f_NO_3)/2;
   f_N2_3 = 1 - f_O2_3 - f_NO_3 - f_O_3 - f_N_3;
   n_N2 = n_N2 * f_N2_3;
   n_O2 = density_f_exc(Tv1, f_O2_3, O2);
   n_NO = density_f_exc(Tv1, f_NO_3, NO);
   y0=[n_N2; n_O2; n_NO; f_N_3; f_O_3; T3];
   options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+1);
   [X, Y]=ode15s(@(t, y) Rpart_0D(t, y, kinetics), xspan, y0, ...
                                                            options_s);
%    [X, Y]=ode15s(@(t, y) Rpart_ODE_0D(t, y, kinetics), xspan, y0, ...
%                                                             options_s);

   t=X*t0;
   Y(:, 1:end-1)=Y(:, 1:end-1)*n0;
   Y(:, end)=Y(:, end)*T0;
   T=Y(:, end);
   Tv = N2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
   out(i_ini, i_vibr, i_U).res=[t, Y, Tv];

   disp('Conservation laws check')
   check_CL_0D([1 1], Y, kinetics, 1);
  end
 end
end

figure
plot(t*1e3, T, t*1e3, Tv, '-.', 'linewidth', 1.5)
legend('T', 'Tv', 'location', 'best')
% xlim([6e-2 1.5e-1])

rmpath('../src/')
toc
end