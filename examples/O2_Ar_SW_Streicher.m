% The main function for the macroparameters calculation behind reflected SW
% for Streicher's experiment conditions in O2-Ar mixture.
% [1] J. Streicher, A. Krish, R. Hanson, Coupled vibration-dissociation 
% time-histories and rate measurements in shock-heated, nondilute O2 and 
% O2-Ar mixtures from 6000 to 14000 K, Physics of Fluids 33 (5) (2021)
%056107. doi:10.1063/5.0048059
% 06.04.2023 Denis Kravchenko

tic
clearvars;
% constants
k=1.380649e-23;             % Boltzmann constant, J/K
Torr=133.322368; 
Na=6.02214076e23;           % Avogadro const
addpath('../src/')
addpath('../data/')
load('particles.mat', 'O2', 'O', 'Ar');
O2.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;
Ar.num_elex_levels=1;

%initialization of structures dat and dat1
tmp.time=0; tmp.T=0; tmp.Tv=0; tmp.nO=0; tmp.nO2=0; tmp.nAr=0; tmp.p=0;
dat(2,4,9)=tmp;

tmp1.time=0; tmp1.T=0; tmp1.Tv=0; tmp1.nO=0; tmp1.nO2=0; tmp1.nAr=0;
tmp1.p=0;
dat1(2,4,9,2)=tmp1;

% initial conditions
init_c=[ %  f;  p0,     Torr;   v0, m/s;    T0, K;   v0_1
        0.5 0.18 2210 296 910
        0.5 0.1 2520 296 1030
        0.5 0.06 2630 296 1080
        0.2 0.47 2000 296 940
        0.2 0.25 2300 296 1070
        0.2 0.06 2670 296 1230
        1 0.13 2220 296 770     % pure O2
        1 0.07 2510 296 870     % pure O2
        1 0.05 2760 296 950     % pure O2
        ];
for i_ini=9 % [1 2 3 4 5 6 7 8 9] % choosing desired initial coonditions
for i_U=4 % [2 3 4]    % choosing desired U dissociation parameter model
% 2 is for D/6k; 3 is for 3T; 4 is for inf
for i_vibr=2 % [1 2]  % choosing vibrational energy exchange model
% 1 is for SSH; 2 is for FHO
for rel= 2 %1:2     % if relaxation between incident and reflected waves 
% frozen? 1 -relaxation off; 2 - relaxation on
    f=init_c(i_ini, 1); %molar fraction of O2
    p0=init_c(i_ini, 2)*Torr; %initial pressure in shock tube
    v0=init_c(i_ini, 3);   % velocity of incident SW, m/s
    v0_r=init_c(i_ini, 5); %velocity of reflected SW, m/s
    T0=init_c(i_ini, 4);  % initial temperature in shock tube, K
    T0buf=T0; %buffer variable for initial temperature
    n0=p0/(k*T0); % initial number density, m-3
    n0buf=n0; %buffer variable for initial number density
    rho0=n0*((1-f)*Ar.mass + f*O2.mass);
    [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0 ,f);
    sigma0 = pi*O2.diameter^2;
    Delta = 1 / sqrt(2) / n0 / sigma0;  %free path length
    % num=0;
    % index{1}=0;
    if f==1
    % Ps={num, O2, O};
        Ps = {O2, O};
    else
    % Ps={num, O2, O, Ar};
        Ps = {O2, O, Ar};
    end
    % for ind=2:length(Ps)
    %     num_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    %     num=num+num_states;
    %     first=index{ind-1}(end)+1;
    %     index{ind}=first:first+num_states-1;
    % end
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
    Reacs_keys={'Diss', 'VT', 'VV'};
    reacs_val={Diss, model_VT, model_VT};
    % kinetics.Ps=Ps(2:end);
    kinetics.Ps = Ps;
    kinetics.num_Ps=length(kinetics.Ps);
    kinetics.index = indexes_for_Ps(Ps);
    kinetics.num_eq = kinetics.index{end}(end);
    % kinetics.num_eq=num;
    kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
    % kinetics.index=index(2:end);
    kinetics.n0=n0;
    kinetics.v0=v0;
    kinetics.T0=T0;
    kinetics.Delta=Delta;
    deltas=5e-3; %distance between transducer and endwall
    timewave=(v0 + v0_r)/(v0*v1)*deltas/v0_r;  %time between SWs in physical 
    % timescale
    x_w=v0*timewave; 
    xspan=[0 x_w]./Delta;
    n=density_f_exc(T0, n1*f, O2); %equilibrium distribution
    y0=zeros(kinetics.num_eq+2, 1); %initial vector
    y0(1:length(n))=n;
    y0(end-1)=v1;
    y0(end)=T1;
    if f~=1
        y0(kinetics.index{end})=n1*(1-f);
    end
    options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-19, ...
                                    'NonNegative', 1:kinetics.num_eq+2); 
    if rel==2  %if relaxation between SWs ON
    [X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), ...
                            xspan, y0, options_s); 
    %incident SW
    X=X*Delta;
    Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
    Y(:, end-1)=Y(:, end-1)*v0;
    Y(:, end)=Y(:, end)*T0;
    n_O=Y(:, length(n)+1);
    if f~=1
        n_Ar=Y(:, kinetics.index{end}); %Argon number density in case with
        % diluted oxygen
    else
        n_Ar=zeros(length(n_O), 1); %Argon number density = 0 in case with
        % pure oxygen
    end
    n_O2=sum(Y(:, 1:length(n)), 2);
    T=Y(:, end);
    p=(n_O+n_O2+n_Ar).*k.*T /Torr; %pressure
    Tv = O2.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
    time_ms=X./v0*1e6;
    elseif rel==1 %if relaxation between SWs is neglected,
    % then initial parameters before reflected SWs are extracted by shock 
    % relations
    Y=y0*n0;
    Y(end)=y0(end)*T0;
    T=Y(end);
    Y(end-1)=v0*y0(end-1);
    Y=Y';
    if f~=1
        n_Ar=Y(kinetics.index{end});
    else
        n_Ar=0;
    end
    n_O=Y(length(n)+1);
    n_O2=sum(Y(1:length(n)));
    end

    rhov0=n0*(f*O2.mass+(1-f)*Ar.mass) * v0;      % rho0*v0
    rhov2p0=n0*(f*O2.mass+(1-f)*Ar.mass)* v0^2 + n0*k*T0; % rho0*v0^2+p0
    e_i=[];
    for ind_e=1:O2.num_elex_levels
        e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
    end
    En0=n0*e_i*n/n1 + n0*k*f*T0 + 1.5*n0*k*T0 + n0*O2.form_e*f;
    Ep0=(En0+n0*k*T0)/(n0*(f*O2.mass+(1-f)*Ar.mass))+0.5*v0^2;      
    % (E0+p0)/rho0+v0^2/2
    if rel==2 %if relaxation ON then checking the conservation laws,
        % otherwise it's pointless
        disp('Conservation laws check, incident SW')
        check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 0);
    end
    
    
    %% REFL
    %calculation initial parameters before reflected SW
    if rel==2 %if relaxation ON then we took last values from previous
    % calculation
    n0=sum(Y(end, 1:end-2),2);   % m-3
    v0=v0+v0_r-Y(end, end-1);   % m/s
    %v0=v0_r+Y(end,end-1);
    T0=Y(end, end);   % K
    rho0_1=n_O2(end)*O2.mass+ n_O(end)*O.mass+n_Ar(end)*Ar.mass;
    elseif rel==1 
    n0=sum(Y(1:end-2));
    v0=v0+v0_r-Y(end-1);
    T0=Y(end);
    rho0_1=n0*(f*O2.mass + (1-f)*Ar.mass);
    end
    [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0_1, n_O2(end)/(n_O2(end) + ...
        n_Ar(end) + n_O(end)));
    Delta = 1 / sqrt(2) / n0 / sigma0;
    kinetics.n0=n0;
    kinetics.v0=v0;
    kinetics.T0=T0;
    kinetics.Delta=Delta;
    timewave=450*1e-6;
    x_w=v0_r*timewave;
    xspan=[0 x_w]./Delta;
    y0_1=zeros(kinetics.num_eq+2, 1);
    if rel==2 
    %if relaxation on then we recalculate distribution from previous
    % calculation
    y0_1(1:end)=Y(end, :).*((1/n0)*n1); 
    elseif rel==1 
    %if relaxation off then we use the equilibrium distribution 
    % with initial! temperature before incident SW
    n=density_f_exc(T0buf, n1*f, O2);
    y0_1(1:length(n))=n;
    if f~=1
        y0_1(kinetics.index{end})=n_Ar/n0*n1;
    end
    end
    y0_1(end-1)=v1;
    y0_1(end)=T1;
        % great for an accurate simulation
    options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-15, ... 
                                    'NonNegative', 1:kinetics.num_eq+2);
        % enough for debugging
    % options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
    %                                 'NonNegative', 1:kinetics.num_eq+2); 
    [X_1, Y_1]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics),...
        xspan, y0_1, options_s);
    X_1=X_1*Delta;
    Y_1(:, 1:end-2)=Y_1(:, 1:end-2)*n0;
    Y_1(:, end-1)=Y_1(:, end-1)*v0;
    Y_1(:, end)=Y_1(:, end)*T0;
    T_1=Y_1(:, end);
    n_O_1=Y_1(:, length(n)+1);
    n_O2_1=sum(Y_1(:, 1:length(n)), 2);
    if f~=1
        n_Ar_1=Y_1(:, kinetics.index{end});
    else
        n_Ar_1=zeros(length(n_O_1),1);
    end
    p_1=(n_O2_1+n_O_1+n_Ar_1)*k.*T_1 /Torr;
    Tv_1 = O2.ev_i{1}(2)./(k*log(Y_1(:,1)./Y_1(:,2)));
    time_ms_1=X_1./v0_r*1e6;
    rhov0_1=rho0_1 * v0;                    % rho0*v0
    rhov2p0_1=rho0_1* v0^2 + n0*k*T0;      % rho0*v0^2+p0
    e_i=[];
    for ind_e=1:O2.num_elex_levels
        e_i=[e_i, O2.ev_i{ind_e}+O2.ev_0(ind_e)+O2.e_E(ind_e)];
    end
    En0_1=n0*e_i*y0_1(1:length(n))/n1 + k*T0*n_O2(end) + 1.5*n0*k*T0 +...
    n_O2(end)*O2.form_e + n_O(end)*O.form_e; 
    %we use the initial distribution here
    Ep0_1=(En0_1+n0*k*T0)/(rho0_1)+0.5*v0^2;    % (E0+p0)/rho0+v0^2/2
    disp('Conservation laws check, reflected SW')
    check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 0);
  
    %This is where the output data is stored. 
    % They contain the evolution of temperatures, number densities,
    % and pressure between the SWs and behind the reflected SW
    if rel==2
        resSt.time=time_ms;
        resSt.T=T;
        resSt.Tv=Tv;
        resSt.nO=n_O;
        resSt.nO2=n_O2;
        if f~=1
            resSt.nAr=n_Ar;
        else
            resSt.nAr=zeros(length(time_ms),1);
        end
        resSt.p=p;
        dat(i_vibr,i_U,i_ini)=resSt;
    end

    resSt_1.time=time_ms_1;
    resSt_1.T=T_1;
    resSt_1.Tv=Tv_1;
    resSt_1.p=p_1;
    resSt_1.nO2=n_O2_1/Na;
    resSt_1.nO=n_O_1/Na; 
    if f~=1
        resSt_1.nAr=n_Ar_1/Na;
    else
        resSt_1.nAr=zeros(length(time_ms_1),1);
    end
    dat1(i_vibr,i_U,i_ini,rel)=resSt_1;
end
end
end
end
%%
%if you want to save your data in .mat file, uncomment following raws
%save(['..\data\O2_Ar Streicher experiment\O2Ar_betweenSWs_output'], 'dat');
%save(['..\data\O2_Ar Streicher experiment\O2Ar_behindRSW_output'], 'dat1');

rmpath('../src/')
rmpath('../data/')
toc                           