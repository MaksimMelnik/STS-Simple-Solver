function out=O2O_SW_Shatalov
% The main function for the macroparameters calculation behind SW for
% Shatalov's experiment conditions.
% 27.12.2022 Maksim Melnik

tic
    % constants
k=1.380649e-23;                     % Boltzmann constant, J/K
addpath('../src/')
addpath('../data/')
load('particles.mat', 'O2', 'O');   % loading particles data
O2.num_elex_levels=1;               % no electronic excitation
O.num_elex_levels=1;

    % initial conditions, should be rewritten using recalculation 
    % functions on SW
init_c=[    % n0, m-3;              v0, m/s;    T0, K
            2.582622078206076e+22   4440        2.991113007914224e+02
            3.116082265979940e+22   4130        3.098808789880962e+02
            3.255838505493098e+22   3950        2.965793020605731e+02
            3.132658068243728e+22	3400        3.082412093964783e+02
            6.950596402478228e+22	3070        2.778507787437617e+02  ];
for i_ini=1 % [1 2 3 4 5] % choosing desired initial coonditions
 for i_U=4 % [2 3 4]      % choosing desired U dissociation parameter model
                          % 2 is for D/6k; 3 is for 3T; 4 is for inf
  for i_vibr=2 % [1 2]    % choosing vibrational energy exchange model
                          % 1 is for SSH; 2 is for FHO
    % 0 --- before SW, 1 --- behind SW
   n0 = init_c(i_ini, 1); % m-3; characteristic number density
   v0 = init_c(i_ini, 2); % m/s; characteristic velocity
   T0 = init_c(i_ini, 3); % K; characteristic temperature
   rho0 = n0 * O2.mass;
        % conservation laws to evaluate macroparameters behind SW
   [n1, v1, T1] = in_con_SW(n0, v0, T0, rho0, 1); % dimentionless numbers
   
   sigma0 = pi*O2.diameter^2;
   Delta = 1 / sqrt(2) / n0 / sigma0;   % characteristic length

num=0;
index{1}=0;
Ps={num, O2, O};          % chemical composition of the mixture
for ind=2:length(Ps)
    num_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_states-1;
end
Diss.Arrhenius='Park';    % parameters in dissociation Arrhenius law
Diss.rec=true;            % is recombination included?
Diss.NEmodel='MT';        % dissociation non-equilibrium model
   switch i_U             % non-equilibrium parameter U in dissociation
    case 2                % model
	 Diss.U='D/6k';
	case 3
	 Diss.U='3T';
	case 4
	 Diss.U='inf';
   end
   switch i_vibr          % vibrational energy exchange model
    case 1
	 model_VT='SSH';
	case 2
	 model_VT='FHO';
   end
Reacs_keys={'Diss', 'VT', 'VV'};
reacs_val={Diss, model_VT, model_VT};
kinetics.Ps=Ps(2:end);    % inicialization of kinetics variable
kinetics.num_Ps=length(kinetics.Ps);
kinetics.num_eq=num;
kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
kinetics.index=index(2:end);
kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;
  xspan=[0 1e0]/Delta;            % interval of integration
  n=density_f_exc(T0, n1, O2);    % vibrational distribution, VDF
  y0=zeros(kinetics.num_eq+2, 1); % initial values behind SW
  y0(1:length(n))=n;
  y0(end-1)=v1;
  y0(end)=T1;
  options_s = odeset('RelTol', 1e-5, ... solving accuracy parameters
                    'AbsTol', 1e-8, 'NonNegative', 1:kinetics.num_eq+2); 
  [X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), ... solving
                                                    xspan, y0, options_s);

X=X*Delta;                      % dimensioning
Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
Y(:, end-1)=Y(:, end-1)*v0;
Y(:, end)=Y(:, end)*T0;
T=Y(:, end);
Tv = O2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2))); % vibrational temperature
   time_ms=X./v0*1e6;                        % time in ms
   out(i_ini, i_vibr, i_U).res=[X, Y, Tv, time_ms]; % output variable

    % checking conservation laws
   rhov0=n0*O2.mass * v0;                    % rho0*v0
   rhov2p0=n0*O2.mass * v0^2 + n0*k*T0;      % rho0*v0^2+p0
   e_i=[];
   for ind_e=1:O2.num_elex_levels
    e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
   end
   En0=n0*e_i*n/n1 + n0*k*T0 + 1.5*n0*k*T0 + n0*O2.form_e;
   Ep0=(En0+n0*k*T0)/(n0*O2.mass)+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
   disp('Conservation laws check')
   check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 0);
  end
 end
end

 % drawing to check if it's fine
figure
semilogx(X, T, X, Tv, 'linewidth', 1.5)

rmpath('../src/')
rmpath('../data/')
toc
end