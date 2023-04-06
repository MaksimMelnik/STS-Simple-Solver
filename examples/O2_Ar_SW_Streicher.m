
% The main function for the macroparameters calculation behind reflected SW for
% Streicher's experiment conditions in O2-Ar mixture.
% 06.04.2023 Denis Kravchenko

tic
    % constants
k=1.380649e-23;  % Boltzmann constant, J/K
Torr=133.322368;
Na=6.02214076e23;
addpath('../src/')
addpath('../data/')
load('O2_O.mat') %#ok<LOAD>
load('particle.mat');
O2.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;

    % initial conditions
init_c=[ % f;  p0, Torr;   v0, m/s;   T0, K;   v0_1
    0.5 0.18 2210 296 910
    0.5 0.1 2520 296 1030
    0.5 0.06 2630 296 1080
    0.2 0.47 2000 296 940
    0.2 0.25 2300 296 1070
    0.2 0.06 2670 296 1230
    ];
for i_ini=[1:6]
 for i_U=[2:4]
  for i_vibr=[1 2]%:2
   f=init_c(i_ini, 1);
   p0=init_c(i_ini, 2)*Torr;
   v0=init_c(i_ini, 3);   % m/s
   v0_r=init_c(i_ini, 5); %m/s
   T0=init_c(i_ini, 4);   % K
   n0=p0/(k*T0);   % m-3
   n0buf=n0;
   NN=in_con_Ar([O2.mass, v0, T0, Ar.mass ,f]);
   n1=NN(1);   % DN
   v1=NN(3);   % DN
   T1=NN(2);   % DN
   
   sigma0 = pi*O2.diameter^2;
Delta = 1 / sqrt(2) / n0 / sigma0;

num=0;
index{1}=0;
Ps={num, O2, O, Ar};
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
deltas=5e-3;
timewave=(v0 + v0_r)/(v0*v1)*deltas/v0_r;
x_w=v0*timewave;
xspan=[0 x_w]./Delta;
n=density_f_exc(T0, n1*f, O2);
y0=zeros(kinetics.num_eq+2, 1);
y0(1:length(n))=n;
y0(end-1)=v1;
y0(end)=T1;
y0(kinetics.index{end})=n1*(1-f);
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
[X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, y0, options_s);

X=X*Delta;
Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
Y(:, end-1)=Y(:, end-1)*v0;
Y(:, end)=Y(:, end)*T0;
n_a=Y(:, length(n)+1);
n_Ar=Y(:, kinetics.index{end});
n_m=sum(Y(:, 1:length(n)), 2);
T=Y(:, end);
Tv = O2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
time_ms=X./v0*1e6;
   out(i_ini, i_vibr, i_U).res=[X, Y, Tv, time_ms];

rhov0=n0*(f*O2.mass+(1-f)*Ar.mass) * v0;                    % rho0*v0
rhov2p0=n0*(f*O2.mass+(1-f)*Ar.mass)* v0^2 + n0*k*T0;      % rho0*v0^2+p0
e_i=[];
for ind_e=1:O2.num_elex_levels
  e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
end
En0=n0*e_i*n/n1 + n0*k*f*T0 + 1.5*n0*k*T0 + n0*O2.form_e*f;
Ep0=(En0+n0*k*T0)/(n0*(f*O2.mass+(1-f)*Ar.mass))+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
disp('Conservation laws check')
check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 1);

%% REFL
n0=sum(Y(end, 1:end-3),2);   % m-3
   v0=v0+v0_r-Y(end, end-1);   % m/s
   T0=Y(end, end);   % K
   NN=in_con_Ar([O2.mass, v0, T0, Ar.mass, f]);
   n1=NN(1);   % DN
   v1=NN(3);   % DN
   T1=NN(2);   % DN
   
Delta = 1 / sqrt(2) / n0 / sigma0;

kinetics.n0=n0;
kinetics.v0=v0;
kinetics.T0=T0;
kinetics.Delta=Delta;
timewave=50*1e-6;
x_w=v0_r*timewave;
xspan=[0 x_w]./Delta;
%n=density_f_exc(T0, f*n1, O2);
y0_1=zeros(kinetics.num_eq+2, 1);
y0_1(1:end)=Y(end, :).*((1/n0)*n1);
y0_1(end-1)=v1;
y0_1(end)=T1;
%y0_1(kinetics.index{end})=Y(end, kinetics.index{end});
 options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
[X_1, Y_1]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, y0_1, options_s);

X_1=X_1*Delta;
Y_1(:, 1:end-2)=Y_1(:, 1:end-2)*n0;
Y_1(:, end-1)=Y_1(:, end-1)*v0;
Y_1(:, end)=Y_1(:, end)*T0;
T_1=Y_1(:, end);
n_a_1=Y_1(:, length(n)+1);
n_Ar_1=Y_1(:, kinetics.index{end});
n_m_1=sum(Y_1(:, 1:length(n)), 2);
p_1=(n_m_1+n_a_1+n_Ar_1)*k.*T_1;
Tv_1 = O2.ev_i{1}(2)./(k*log(Y_1(:,1)./Y_1(:,2)));
time_ms_1=X_1./v0_r*1e6;
   out(i_ini, i_vibr, i_U).res=[X_1, Y_1, Tv_1, time_ms_1];

rhov0_1=n0*(f*O2.mass+(1-f)*Ar.mass) * v0;                    % rho0*v0
rhov2p0_1=n0*(f*O2.mass+(1-f)*Ar.mass)* v0^2 + n0*k*T0;      % rho0*v0^2+p0
e_i=[];
for ind_e=1:O2.num_elex_levels
  e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
end
En0_1=n0*e_i*n/n1 + n0*k*T0*f + 1.5*n0*k*T0 + n0*O2.form_e*f;
Ep0_1=(En0_1+n0*k*T0)/(n0*(f*O2.mass+(1-f)*Ar.mass))+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
disp('Conservation laws check')
check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 1);


resSt.time=time_ms;
resSt.T=T;
resSt.Tv=Tv;
resSt.na=n_a;
resSt.nm=n_m;
resSt.nAr=n_Ar;
dat(i_vibr,i_U,i_ini)=resSt;


resSt_1.time=time_ms_1;
resSt_1.T=T_1;
resSt_1.Tv=Tv_1;
resSt_1.p=p_1/Torr;
resSt_1.nm_n=n_m_1/Na;
resSt_1.na=n_a_1/Na; 
dat1(i_vibr,i_U,i_ini)=resSt_1;
  end
 end
end

%%
% figure("Position", [0, 0, 700, 330])
% tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight")
% nexttile
% plot(time_ms, T, time_ms, Tv, 'linewidth', 1.5);
% xlim([0 50]);
% nexttile
% plot(time_ms_1, T_1, time_ms_1, Tv_1, 'linewidth', 1.5);
% xlim([0 50]);
% 
% load('dat.mat');
% load('data_new.mat');
% figure("Position", [0, 0, 700, 330])
% tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight")
% nexttile
% plot(dat(1,3,3,1,2).time, dat(1,3,3,1,2).T, dat(1,3,3,1,2).time, dat(1,3,3,1,2).Tv);
% nexttile
% plot(data_1(1,3,3,1,2).time, data_1(1,3,3,1,2).T, data_1(1,3,3,1,2).time, data_1(1,3,3,1,2).Tv01);
% xlim([0 50]);

save(['datAr_02.mat'], 'dat');
save(['datAr_1_02.mat'], 'dat1');

rmpath('../src/')
rmpath('../data/')
toc
