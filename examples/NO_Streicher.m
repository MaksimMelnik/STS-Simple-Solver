% The main function for the macroparameters calculation behind reflected SW
% for Streicher's experiment conditions in NO-Ar mixture.
% 06.04.2023 Denis Kravchenko
tic
clearvars;
% constants
k=1.380649e-23;  % Boltzmann constant, J/K
Torr=133.322368;
Na=6.02214076e23;
addpath('../src/')
addpath('../data/')
load('particles.mat', "NO", "N", "O", "Ar", "O2", "N2");
NO.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;
N.num_elex_levels=1;
O2.num_elex_levels=1;
N2.num_elex_levels=1;
%O2.num_vibr_levels=1;
%N2.num_vibr_levels=1;

% initial conditions
init_c=[ % f;  p0, Torr;   v0, m/s;   T0, K;   v0_1
    0.02 4.31 1265 296 668
    0.02 1.4 1578 296 814
    0.02 0.36 1801 296 921
    0.02 0.35 2004 296 1019
    0.01 10.83 1015 296 557
    0.01 4.71 1244 296 661
    0.01 1.08 1591 296 824
    0.01 0.89 1797 296 923
    0.004 14.95 973 296 539
    0.004 5.35 1241 296 662
    0.004 2.34 1597 296 830
    ];
for i_ini=1:11
for i_U=2:4 %choosing desired U dissociation parameter model
%2 is for D/6k; 3 is for 3T; 4 is for inf
for i_vibr=2% choosing desired vibrational energy exchange model
%1 for SSH; 2 for FHO
for rel=2 % 1 -relaxation off; 2 - relaxation on
    f=init_c(i_ini, 1); %molar fraction of NO
    p0=init_c(i_ini, 2)*Torr; %initial pressure in shock tube
    v0=init_c(i_ini, 3);   % velocity of incident SW, m/s
    v0_i=v0;
    v0_r=init_c(i_ini, 5); %velocity of reflected SW, m/s
    T0=init_c(i_ini, 4); % initial temperature in shock tube, K
    T0buf=T0; %buffer variable for initial temperature
    n0=p0/(k*T0);   % initial number density, m-3
    n0buf=n0; %buffer variable for initial number density
%     NN=in_con_Ar([NO.mass, v0, T0, Ar.mass ,f]); %dimensionless variables
%     n1=NN(1);   % DN
%     v1=NN(3);   % DN
%     T1=NN(2);   % DN

    rho0=n0*((1-f)*Ar.mass + f*NO.mass);
    [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0 ,f);

    sigma0 = pi*NO.diameter^2;
    Delta = 1 / sqrt(2) / n0 / sigma0; %free path length
    num=0;
    index{1}=0;
    Ps={num, NO, O, N, O2, N2, Ar};
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
    Exch=1;
    Reacs_keys={'VT','VV'}; %chemical relaxation and dissociation between 
    %SWs is negligible, only VT and VV processes 
    %Reacs_keys={'Diss', 'VT', 'VV'};
    reacs_val={model_VT, model_VT};
    %reacs_val={Diss,  model_VT, model_VT};
    kinetics.Ps=Ps(2:end);
    kinetics.num_Ps=length(kinetics.Ps);
    kinetics.num_eq=num;
    kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
    kinetics.index=index(2:end);
    kinetics.n0=n0;
    kinetics.v0=v0;
    kinetics.T0=T0;
    kinetics.Delta=Delta;
    
    %determine index numbers of molecules
    names=repmat("", length(kinetics.Ps), 1);
    serial_index=zeros(length(kinetics.Ps), 1);
    for i=1:length(kinetics.Ps)
        names(i)=string(kinetics.Ps{i}.name);
        serial_index(i)=i;
        IndexOfMolecules=containers.Map(names,serial_index);
    end
    kinetics.IndexOfMolecules=IndexOfMolecules;
    deltas=5e-3; %distance from the transducer to the endwall
    timewave=(v0 + v0_r)/(v0*v1)*deltas/v0_r; %time between passing the SWs
    x_w=v0*timewave;
    xspan=[0 x_w]./Delta;
    n=density_f_exc(T0, n1*f, NO);
    y0=zeros(kinetics.num_eq+2, 1);
    y0(1:length(n))=n;
    y0(end-1)=v1;
    y0(end)=T1;
    y0(kinetics.index{end})=n1*(1-f);
    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
    'NonNegative', 1:kinetics.num_eq+2);
    if rel==2 %if relaxation between SWs on
    [X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, ...
        y0, options_s);
    X=X*Delta;
    Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
    Y(:, end-1)= Y(:, end-1)*v0 ;
    Y(:, end)=Y(:, end)*T0;
    n_O=Y(:, kinetics.index{IndexOfMolecules("O")});
    n_N=Y(:, kinetics.index{IndexOfMolecules("N")});
    n_O2 = sum(Y(:, kinetics.index{IndexOfMolecules("O2")}), 2);
    n_N2 = sum(Y(:,kinetics.index{IndexOfMolecules("N2")}),2);
    n_Ar=Y(:, kinetics.index{IndexOfMolecules("Ar")});
    n_NO=sum(Y(:, kinetics.index{IndexOfMolecules("NO")}), 2);
    T=Y(:, end);
    p=(n_O+n_N + n_NO +n_Ar + n_O2 + n_N2).*k.*T /Torr;
    Tv = NO.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
    time_ms=X./v0*1e6;
    elseif rel==1 
    %if relaxation between SWs is off, than using the R-H relation
    Y=y0*n0;
    Y(end)=y0(end)*T0;
    T=Y(end);
    Y(end-1)=v0*y0(end-1);
    Y=Y';
    n_Ar=Y(kinetics.index{IndexOfMolecules("Ar")});
    n_O=Y(kinetics.index{IndexOfMolecules("O")});
    n_N=Y(kinetics.index{IndexOfMolecules("N")});
    n_NO=sum(Y(kinetics.index{IndexOfMolecules("NO")}));
    end
    rhov0=n0*(f*NO.mass+(1-f)*Ar.mass) * v0;                   
    % rho0*v0
    rhov2p0=n0*(f*NO.mass+(1-f)*Ar.mass)* v0^2 + n0*k*T0;      
    % rho0*v0^2+p0
    e_i=[];
    for ind_e=1:NO.num_elex_levels
    e_i=[e_i, NO.ev_i{ind_e}+NO.ev_0(ind_e)+NO.e_E(ind_e)];
    end
    En0=n0*e_i*n/n1 + n0*k*f*T0 + 1.5*n0*k*T0 + n0*NO.form_e*f;
    Ep0=(En0+n0*k*T0)/(n0*(f*NO.mass+(1-f)*Ar.mass))+0.5*v0^2;  
    % (E0+p0)/rho0+v0^2/2
    if rel==2
    disp('Conservation laws check behind ISW')
    check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 0);
    end
    
    %% REFL
   % taking into account exchange reactions
    Reacs_keys={'Diss','Exch', 'VT', 'VV'};
    reacs_val={Diss, Exch,  model_VT, model_VT};

    %without taking into account exchange reactions
    %Reacs_keys={'Diss', 'VT', 'VV'};
    %reacs_val={Diss,  model_VT, model_VT};
    kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
    if rel==2
    n0=sum(Y(end, 1:end-2),2);   % m-3
    v0=v0+v0_r-Y(end, end-1);   % m/s
    T0=Y(end, end);   % K
    elseif rel==1
    n0=sum(Y(1:end-2));
    v0=v0+v0_r-Y(end-1);
    T0=Y(end);
    end
    NN=in_con_Ar([NO.mass, v0, T0, Ar.mass, f]);
    n1=NN(1);   % DN
    v1=NN(3);   % DN
    T1=NN(2);   % DN
    
    kinetics.n0=n0;
    kinetics.v0=v0;
    kinetics.T0=T0;
    kinetics.Delta=Delta;
    timewave=1000*1e-6;
    x_w=v0_r*timewave;
    xspan=[0 x_w]./Delta;
    %n=density_f_exc(T0, f*n1, O2);
    y0_1=zeros(kinetics.num_eq+2, 1);
    if rel==2
    y0_1(1:end)=Y(end, :).*((1/n0)*n1);
    elseif rel==1
    n=density_f_exc(T0buf, n1*f, NO);
    y0_1(1:length(n))=n;
    y0_1(kinetics.index{end})=n_Ar/n0*n1;
    end
    y0_1(end-1)=v1;
    y0_1(end)=T1;
    %y0_1(kinetics.index{end})=Y(end, kinetics.index{end});
    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
    'NonNegative', 1:kinetics.num_eq+2);
    [X_1, Y_1]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics),...
        xspan, y0_1, options_s);
    
    
    X_1=X_1*Delta;
    Y_1(:, 1:end-2)=Y_1(:, 1:end-2)*n0;
    Y_1(:, end-1)=Y_1(:, end-1)*v0;
    Y_1(:, end)=Y_1(:, end)*T0;
    T_1=Y_1(:, end);
    n_O_1=Y_1(:,  kinetics.index{IndexOfMolecules("O")});
    n_N_1=Y_1(:,  kinetics.index{IndexOfMolecules("N")});
    n_Ar_1=Y_1(:, kinetics.index{IndexOfMolecules("Ar")});
    n_NO_1=sum(Y_1(:, kinetics.index{IndexOfMolecules("NO")}), 2);
    n_O2_1 = sum(Y_1(:, kinetics.index{IndexOfMolecules("O2")}), 2);
    n_N2_1 = sum(Y_1(:, kinetics.index{IndexOfMolecules("N2")}),2);
    p_1=(n_NO_1+n_O_1+ + n_N_1 + n_Ar_1 + n_O2_1 + n_N2_1)*k.*T_1 /Torr;
    Tv_NO_1 = NO.ev_i{1}(2)./(k*log(Y_1(:, ...
        kinetics.index{IndexOfMolecules("NO")}(1))./...
        Y_1(:,kinetics.index{IndexOfMolecules("NO")}(2))));
    Tv_O2=O2.ev_i{1}(2)./(k*log(Y_1(:, ...
        kinetics.index{IndexOfMolecules("O2")}(1))./...
     Y_1(:,kinetics.index{IndexOfMolecules("O2")}(2))));
    Tv_N2=N2.ev_i{1}(2)./(k*log(Y_1(:, ...
        kinetics.index{IndexOfMolecules("N2")}(1))./...
     Y_1(:,kinetics.index{IndexOfMolecules("N2")}(2))));
    time_ms_1=X_1./v0_r*1e6;
    if rel==2
    rho0=n_NO(end)*NO.mass + n_O(end)*O.mass + n_N(end)*N.mass + ...
        n_O2(end)*O2.mass + n_N2(end)*N2.mass + n_Ar(end)*Ar.mass;
    else
    rho0=n0*NO.mass*f + n0*Ar.mass*(1-f);
    end
    rhov0_1=rho0 * v0;                    % rho0*v0
    rhov2p0_1=rho0* v0^2 + n0*k*T0;      % rho0*v0^2+p0
    e_i=[];
    for ind_e=1:NO.num_elex_levels
    e_i=[e_i, NO.ev_i{ind_e}+NO.ev_0(ind_e)+NO.e_E(ind_e)];
    end
    En0_1=n0*e_i*y0_1(1:length(n))/n1 + k*T0*n_NO(end) + 1.5*n0*k*T0 +...
    n_NO(end)*NO.form_e + n_N(end)*N.form_e+n_O(end)*O.form_e;
    Ep0_1=(En0_1+n0*k*T0)/rho0+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
    disp('Conservation laws check behind RSW')
    check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 0);

    %This is where the output data is stored. 
    % They contain the evolution of temperatures, number densities,
    % and pressure between the SWs and behind the reflected SW
    if rel==2
    resSt.time=time_ms;
    resSt.T=T;
    resSt.Tv=Tv;
    resSt.nO=n_O;
    resSt.nN=n_N;
    resSt.nNO=n_NO;
    resSt.nAr=n_Ar;
    resSt.nO2=n_O2;
    resSt.nN2=n_N2;
    resSt.p=p;
    dat(i_vibr,i_U,i_ini)=resSt;
    end
    resSt_1.time=time_ms_1;
    resSt_1.T=T_1;
    resSt_1.TvNO=Tv_NO_1;
    resSt_1.TvO2=Tv_O2;
    resSt_1.TvN2=Tv_N2;
    resSt_1.ni_NO=Y_1(:, kinetics.index{IndexOfMolecules("NO")});
    resSt_1.p=p_1;
    resSt_1.nNO=n_NO_1/Na;
    resSt_1.nO2=n_O2_1/Na;
    resSt_1.nN2=n_N2_1/Na;
    resSt_1.nO=n_O_1/Na;
    resSt_1.nN=n_O_1/Na;
    resSt_1.nAr=n_Ar_1/Na;
    dat1(i_vibr,i_U,i_ini,rel)=resSt_1;
end
end
end

end

%%
%if you want to save your data in .mat file, uncomment following raws
save(['NO_between_SWs_withexch_VDOP1_PARK.mat'], 'dat');
save(['NO_behind_ReflSW_withexch_VDOP1_PARK.mat'], 'dat1');
rmpath('../src/')
rmpath('../data/')
toc                