
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
init_c=[ % n0, m-3;   v0, m/s;   T0, K;   v0_1
    4.240873665727726e+21 2220 296 770
    ];
for i_ini=1
 for i_U=3
  for i_vibr=1%:2
   n0=init_c(i_ini, 1);   % m-3
   n0buf=n0;
   v0=init_c(i_ini, 2);   % m/s
   v0_r=init_c(i_ini, 4); %m/s
   T0=init_c(i_ini, 3);   % K
   NN=in_con_O2([O2.mass, v0, T0]);
   n1=NN(1);   % DN
   v1=NN(3);   % DN
   T1=NN(2);   % DN
   
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
kinetics.num_eq=num + 3;
kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
kinetics.index=index(2:end);
kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;
deltas=5e-3;
timewave=(v0 + v0_r)/(v0*v1)*deltas/v0_r;
x_w=v0*timewave;
xspan=[0 x_w]./Delta;
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

%% REFL
n0=sum(Y(end, 1:end-2),2);   % m-3
   v0=v0+v0_r-Y(end, end-1);   % m/s
   T0=Y(end, end);   % K
   NN=in_con_O2([O2.mass, v0, T0]);
   n1=NN(1);   % DN
   v1=NN(3);   % DN
   T1=NN(2);   % DN
   
Delta = 1 / sqrt(2) / n0 / sigma0;

kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;
timewave=50*1e-6;
x_w=v0*timewave;
xspan=[0 x_w]./Delta;
n=density_f_exc(T0, n1, O2);
y0_1=zeros(kinetics.num_eq+2, 1);
y0_1(1:length(n))=Y(end, 1:length(n)).*((1/n0)*n1);
y0_1(end-1)=v1;
y0_1(end)=T1;
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
[X_1, Y_1]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, y0_1, options_s);

X_1=X_1*Delta;
Y_1(:, 1:end-2)=Y_1(:, 1:end-2)*n0;
Y_1(:, end-1)=Y_1(:, end-1)*v0;
Y_1(:, end)=Y_1(:, end)*T0;
T_1=Y_1(:, end);
Tv_1 = O2.ev_i{1}(2)./(k*log(Y_1(:,1)./Y_1(:,2)));
time_ms_1=X_1./v0_r*1e6;
   out(i_ini, i_vibr, i_U).res=[X_1, Y_1, Tv_1, time_ms_1];

rhov0_1=n0*O2.mass * v0;                    % rho0*v0
rhov2p0_1=n0*O2.mass * v0^2 + n0*k*T0;      % rho0*v0^2+p0
e_i=[];
for ind_e=1:O2.num_elex_levels
  e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
end
En0_1=n0*e_i*n/n1 + n0*k*T0 + 1.5*n0*k*T0 + n0*O2.form_e;
Ep0_1=(En0_1+n0*k*T0)/(n0*O2.mass)+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
disp('Conservation laws check')
check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 0);

  end
 end
end

%%
figure("Position", [0, 0, 700, 330])
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight")
nexttile
plot(time_ms, T, time_ms, Tv, 'linewidth', 1.5);
xlim([0 50]);
nexttile
plot(time_ms_1, T_1, time_ms_1, Tv_1, 'linewidth', 1.5);
xlim([0 50]);

load('dat.mat');
load('data_new.mat');
figure("Position", [0, 0, 700, 330])
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight")
nexttile
plot(dat(1,3,3,1,2).time, dat(1,3,3,1,2).T, dat(1,3,3,1,2).time, dat(1,3,3,1,2).Tv);
nexttile
plot(data_1(1,3,3,1,2).time, data_1(1,3,3,1,2).T, data_1(1,3,3,1,2).time, data_1(1,3,3,1,2).Tv01);
xlim([0 50]);
rmpath('../src/')
rmpath('../data/')
toc
