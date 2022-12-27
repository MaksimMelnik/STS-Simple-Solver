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

    % initial variables just for one case for tests
n0=2.582622078206076e+22;   % m-3
v0=4440;                    % m/s
T0=2.991113007914224e+02;   % K
n1=5.838863656521909e+00;   % DN
v1=1.712662015806821e-01;   % DN
T1=3.617268254080651e+01;   % DN

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
Diss.U='D/6k';
Diss.U='3T';
% Diss.U='inf';
model_VT='FHO';
model_VT='SSH';
Reacs_keys={'Diss', 'VT'};  % VV?
reacs_val={Diss, model_VT};
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
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
[X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, y0, options_s);

X=X*Delta;
Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
Y(:, end-1)=Y(:, end-1)*v0;
Y(:, end)=Y(:, end)*T0;
T=Y(:, end);
Tv = O2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
time_ms=X./v0*1e6;
out=[X, Y, Tv, time_ms];

rhov0=n0*O2.mass * v0;                    % rho0*v0
rhov2p0=n0*O2.mass * v0^2 + n0*k*T0;      % rho0*v0^2+p0
e_i=[];
for ind_e=1:O2.num_elex_levels
  e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
end
En0=n0*e_i*n/n1 + n0*k*T0 + 1.5*n0*k*T0 + n0*O2.form_e;
Ep0=(En0+n0*k*T0)/(n0*O2.mass)+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
disp('Conservation laws check')
check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics);

figure
semilogx(X, T, X, Tv, 'linewidth', 1.5)

rmpath('../src/')
rmpath('../data/')
toc
end