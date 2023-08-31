function out = O2O_SW_Shatalov
% The main function for the macroparameters calculation behind SW for
% Shatalov's experiment conditions [1].
% [1] Ibraguimova et al J. Chem. Phys. 139, 034317 (2013).
% 27.12.2022 Maksim Melnik
tic                                 % calculation time measuring
    % constants
k = 1.380649e-23;                   % Boltzmann constant, J/K
Torr = 133.322368;                  % how much Pa in a one Torr
addpath('../src/')                  % adding the path to the source files
load('../data/particles.mat', 'O2', 'O'); % loading particles data
O2.num_elex_levels=1;               % no electronic excitation
O.num_elex_levels=1;
    % initial condition from Shatalov's paper [1]
init_c_Shatalov = [ % p0, Torr;     T1, K;      v0, m/s;
                    0.8             10820       4440
                    1               9410        4130
                    1               8620        3950
                    1               6470        3400
                    2               5300        3070];
	% initialization of the gas mixture
Ps    = {O2, O};                    % chemical composition of the mixture
index = cell(1, length(Ps));
first = 1;
for ind = 1:length(Ps)
    num_states = sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    index{ind} = first : (first+num_states-1);
    first = index{ind}(end) + 1;
end

for i_ini = 1 % [1 2 3 4 5]         % choosing desired initial coonditions
 cond2 = par_shatalov_f(init_c_Shatalov(i_ini, :));
  % 0 --- before SW, 1 --- behind SW
 v0 = init_c_Shatalov(i_ini, 3);      % m/s; characteristic velocity
 T0 = cond2(2);                       % K; characteristic temperature
 p0 = init_c_Shatalov(i_ini, 1)*Torr; % m-3; characteristic number density
 n0 = p0/k/T0;
 rho0 = n0 * O2.mass;
 sigma0 = pi*O2.diameter^2;
 Delta = 1 / sqrt(2) / n0 / sigma0;   % characteristic length
  % conservation laws to evaluate macroparameters behind SW
 [n1, v1, T1] = in_con_SW(n0, v0, T0, rho0, 1); % dimentionless numbers
 
 for i_U = 4 % [2 3 4]   % choosing desired U dissociation parameter model
                          %  2 is for D/6k; 3 is for 3T; 4 is for inf
  for i_vibr = 2 % [1 2 3]  % choosing vibrational energy exchange model
                          %  1 is for SSH; 2 is for FHO
                          
   Diss.Arrhenius = 'Park';% parameters in dissociation Arrhenius law
   Diss.rec = true;       % is recombination included?
   Diss.NEmodel = 'MT';   % dissociation non-equilibrium model
   switch i_U             % non-equilibrium parameter U in
    case 2                %     the dissociation model
	 Diss.U = 'D/6k';
	case 3
	 Diss.U = '3T';
	case 4
	 Diss.U = 'inf';
   end
   switch i_vibr          % vibrational energy exchange model
    case 1
	 model_VT = 'SSH';
	case 2
	 model_VT = 'FHO';
	case 3
	 model_VT = 'FHO-FR';
   end
   Reacs_keys = {'Diss', 'VT', 'VV'};
   model_VV = model_VT;
   if model_VV == "FHO-FR"
       model_VV = "FHO";
   end
   reacs_val = {Diss, model_VT, model_VV};
   kinetics.Ps = Ps;                 % inicialization of kinetics variable
   kinetics.num_Ps = length(kinetics.Ps);
   kinetics.index = index;
   kinetics.num_eq = kinetics.index{end}(end);
   kinetics.reactions = containers.Map(Reacs_keys, reacs_val);
   kinetics.n0 = n0;
   kinetics.v0 = v0;
   kinetics.T0 = T0;
   kinetics.Delta = Delta;
   xspan = [0 1e0] / Delta;         % interval of integration
   n = density_f_exc(T0, n1, O2);   % initial vibrational distribution,VDF
   y0 = zeros(kinetics.num_eq + 2, 1); % initial values behind SW
   y0(1:length(n)) = n;
   y0(end-1) = v1;
   y0(end) = T1;
   options_s = odeset('RelTol', 1e-5, ... solving accuracy parameters
                    'AbsTol', 1e-8, 'NonNegative', 1:kinetics.num_eq+2); 
   [X, Y] = ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), ... solving
                                                    xspan, y0, options_s);
                                                
   X = X * Delta;                   % dimensioning
   Y(:, 1:end-2) = Y(:, 1:end-2) * n0;
   Y(:, end-1)  =  Y(:, end-1)  *  v0;
   Y(:, end)   =   Y(:, end)   *   T0;
   T = Y(:, end);
        % vibrational temperature
   Tv = O2.ev_i{1}(2) ./ (k * log(Y(:,1)./Y(:,2))); 
   time_mus = X ./ v0 * 1e6;         % time in microseconds
   out(i_ini, i_vibr, i_U).res = [X, Y, Tv, time_mus]; % output variable
   out(i_ini, i_vibr, i_U).time_mus = time_mus;
   out(i_ini, i_vibr, i_U).Tv = Tv;
   out(i_ini, i_vibr, i_U).T  = T;

        % checking conservation laws
   rhov0  =  rho0 * v0;                        % rho0*v0
   rhov2p0 = rho0 * v0^2 + n0*k*T0;          % rho0*v0^2+p0
   e_i = [];
   for ind_e = 1:O2.num_elex_levels
    e_i = [e_i, O2.ev_i{ind_e} + O2.ev_0(ind_e) + O2.e_E(ind_e)];
   end
   En0 = n0*e_i*n/n1 + n0*k*T0 + 1.5*n0*k*T0 + n0*O2.form_e;
   Ep0 = (En0 + n0*k*T0)/rho0 + 0.5*v0^2;      % (E0+p0)/rho0+v0^2/2
   disp('Conservation laws check')
   check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 0);
%    check_CL_SW 0th parameters may be switched to [rho0, v0, T0, En0]
  end
 end
end

 % drawing to check if it's fine
figure
semilogx(X, T, X, Tv, 'linewidth', 1.5)

rmpath('../src/')
toc
end