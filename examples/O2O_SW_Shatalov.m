function out=O2O_SW_Shatalov
% The main function for the macroparameters calculation behind SW for
% Shatalov's experiment conditions.
% 27.12.2022 Maksim Melnik

tic
    % constants
k=1.380649e-23;             % Boltzmann constant, J/K
addpath('../src/')
addpath('../data/')
load('O2_O.mat') %#ok<LOAD>
O2.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;

    % initial conditions
init_c=[ % n0, m-3;   v0, m/s;   T0, K;   n1, DN;   v1, DN;   T1, DN
    2.582622078206076e+22   4440    2.991113007914224e+02 ...
     5.838863656521909e+00   1.712662015806821e-01   3.617268254080651e+01
    3.116082265979940e+22   4130    3.098808789880962e+02 ...
     5.808077724421962e+00  1.721740044550666e-01   3.036555281554731e+01
    3.255838505493098e+22 3.950000000000000e+03 2.965793020605731e+02 ...
     5.799490578609251e+00  1.724289377567720e-01   2.906382663628844e+01
    3.132658068243728e+22	3400    3.082412093964783e+02 ...
     5.722467108795348e+00  1.747498030111898e-01   2.098940479344222e+01
    6.950596402478228e+22	3070	2.778507787437617e+02 ...
     5.694648485891184e+00  1.756034639324195e-01   1.907439874030552e+01
     ];
for i_ini=1:5
 for i_U=2:4
  for i_vibr=2%:2
   n0=init_c(i_ini, 1);   % m-3
   v0=init_c(i_ini, 2);   % m/s
   T0=init_c(i_ini, 3);   % K
   n1=init_c(i_ini, 4);   % DN
   v1=init_c(i_ini, 5);   % DN
   T1=init_c(i_ini, 6);   % DN
   
   sigma0 = pi*O2.diameter^2;
Delta = 1 / sqrt(2) / n0 / sigma0;

num=0;
index{1}=0;
Ps={num, O2, O};
for ind=2:length(Ps)
    num_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_states-1;
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
Reacs_keys={'Diss', 'VT'};  % VV?
Reacs_keys={'Diss', 'VT', 'VV'};
reacs_val={Diss, model_VT};
reacs_val={Diss, model_VT, model_VT};
kinetics.Ps=Ps(2:end);
kinetics.num_Ps=length(kinetics.Ps);
kinetics.num_eq=num;
kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
kinetics.index=index(2:end);
kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;
xspan=[0 1e0]/Delta;
n=density_f_exc(T0, n1, O2);
y0=zeros(kinetics.num_eq+2, 1);
y0(1:length(n))=n;
y0(end-1)=v1;
y0(end)=T1;
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
[X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, y0, options_s);

X=X*Delta;
Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
Y(:, end-1)=Y(:, end-1)*v0;
Y(:, end)=Y(:, end)*T0;
T=Y(:, end);
Tv = O2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
time_ms=X./v0*1e6;
   out(i_ini, i_vibr, i_U).res=[X, Y, Tv, time_ms];

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

figure
semilogx(X, T, X, Tv, 'linewidth', 1.5)

rmpath('../src/')
rmpath('../data/')
toc
end