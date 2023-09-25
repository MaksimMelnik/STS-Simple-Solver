% The main function for the macroparameters calculation behind reflected SW
% for Streicher's experiment conditions in NO-Ar mixture.
% 06.04.2023 Denis Kravchenko

tic
clearvars;
% constants
k=1.380649e-23;  % Boltzmann constant, J/K
Torr=133.322368;
Na=6.02214076e23;

%initialization of structures dat and dat1
tmp.time=0; tmp.T=0; tmp.Tv=0; tmp.nO=0; tmp.nN=0; tmp.nNO=0; tmp.nAr=0;
tmp.nO2=0; tmp.nN2=0; tmp.p=0;
dat(2,4,12,3)=tmp;

tmp1.time=0; tmp1.T=0; tmp1.TvNO=0; tmp1.TvO2=0; tmp1.TvN2=0; tmp1.ni_NO=0;
tmp1.ni_O2=0; tmp1.ni_N2=0; tmp1.p=0; tmp1.nNO=0; tmp1.nO2=0; tmp1.nN2=0;
tmp1.nO=0; tmp1.nN=0; tmp1.nAr=0;
dat1(2,4,12,2,3)=tmp1;

clear tmp tmp1;
addpath('../src/')
addpath('../data/')
load('particles.mat', "NO", "N", "O", "Ar", "O2", "N2");
load('../data/reactions.mat');
ReactZel_1 = Reactions("N2 + O -> NO + N");
ReactZel_2 = Reactions("O2 + N -> NO + O");
Exch = [ReactZel_1("Kunova"), ReactZel_2("Kunova")];

NO.num_elex_levels=1;       % no electronic excitation
O.num_elex_levels=1;
N.num_elex_levels=1;
O2.num_elex_levels=1;
N2.num_elex_levels=1;

% initial conditions
init_c=[ % f;  p0, Torr;   v0, m/s;   T0, K;   v0_1
    0.02 1.52 1768 296 636    % 2% NO; 98% N2
    0.02 0.15 2337 296 815    % 2% NO; 98% N2
    0.02 0.09 2509 296 869    % 2% NO; 98% N2
    0.004 6.09 1243 296 481   % 0.4% NO; 99.6% N2
    0.004 3.40 1767 296 636   % 0.4% NO; 99.6% N2
    0.004 0.59 2153 296 756   % 0.4% NO; 99.6% N2
    0.02 2.06 1526 296 651    % 2% NO; 49% N2; 49% Ar
    0.02 1.2 1995 296 828     % 2% NO; 49% N2; 49% Ar
    0.02 0.22 2297 296 944    % 2% NO; 49% N2; 49% Ar
    0.004 9.68 1064 296 487   % 0.4% NO; 49.8% N2; 49.8% Ar
    0.004 4.63 1470 296 633   % 0.4% NO; 49.8% N2; 49.8% Ar
    0.004 1.18 1959 296 817   % 0.4% NO; 49.8% N2; 49.8% Ar
    ];
for i_exch=2 % [1 2 3]
% 1 - exchange reactions off; 2 - exchange reactions on
% 3 - exchange reaction with average Kunova model and disabled NO vibr.
% spectrum

for i_ini=9 % [1 2 3 4 5 6 7 8 9 10 11 12]
%choosing testcase

for i_U=2 % [2 3 4]
%choosing desired U dissociation parameter model
%2 is for D/6k; 3 is for 3T; 4 is for inf

for i_vibr=2 %[1 2]
% choosing desired vibrational energy exchange model 1 for SSH; 2 for FHO

for i_rel=2 %[1 2]
% 1 -relaxation off; 2 - relaxation on


    if i_exch==3
    NO.num_vibr_levels=1;
    NO.ev_0(1) = 0;
    NO.ev_i{1}=0;
    Exch = [ReactZel_1("Kunova, NO avg"), ReactZel_2("Kunova, NO avg")];
    end
    f=init_c(i_ini, 1); %molar fraction of NO
    p0=init_c(i_ini, 2)*Torr; %initial pressure in shock tube
    v0=init_c(i_ini, 3);   % velocity of incident SW, m/s
    v0_i=v0;
    v0_r=init_c(i_ini, 5); %velocity of reflected SW, m/s
    T0=init_c(i_ini, 4); % initial temperature in shock tube, K
    T0buf=T0; %buffer variable for initial temperature
    n0=p0/(k*T0);   % initial number density, m-3
    n0buf=n0; %buffer variable for initial number density

    if (i_ini<=6)
        rho0=n0*((1-f)*N2.mass + f*NO.mass);
        [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0 ,1); 
        %dimensionless variables
    else
        rho0=n0*((1-f)/2*N2.mass + (1-f)/2*Ar.mass + f*NO.mass);
        [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0 ,f + (1-f)/2);
        %dimensionless variables
    end
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

    Reacs_keys={'VT','VV'}; %dissociation between 
    reacs_val={model_VT, model_VT};
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
    y0=zeros(kinetics.num_eq+2, 1);
    y0(end-1)=v1;
    y0(end)=T1;

    n=density_f_exc(T0, n1*f, NO);
    y0(kinetics.index{IndexOfMolecules("NO")})=n;

    if (i_ini<=6)
        n_=density_f_exc(T0, n1*(1-f), N2);
        y0(kinetics.index{IndexOfMolecules("Ar")})=0;
    else
        n_=density_f_exc(T0, n1*(1-f)/2, N2);
        y0(kinetics.index{IndexOfMolecules("Ar")})=n1*(1-f)/2;
    end
    y0(kinetics.index{IndexOfMolecules("N2")})=n_;
    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
    'NonNegative', 1:kinetics.num_eq+2);
    if i_rel==2 %if relaxation between SWs on
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
    if i_exch~=3 %if i_exch==3 then TvNO has no meaning
    Tv = NO.ev_i{1}(2)./(k*log(Y(:,1)./Y(:,2)));
    end
    time_ms=X./v0*1e6;
    elseif i_rel==1 
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
    n_O2 = sum(Y(:, kinetics.index{IndexOfMolecules("O2")}));
    n_N2 = sum(Y(:,kinetics.index{IndexOfMolecules("N2")}));
    end

        rhov0=rho0 * v0;                   
        % rho0*v0
        rhov2p0=rho0* v0^2 + n0*k*T0;      
        % rho0*v0^2+p0
        e_i=[];
        e_i_=[];
        for ind_e=1:NO.num_elex_levels
        e_i=[e_i, NO.ev_i{ind_e}+NO.ev_0(ind_e)+NO.e_E(ind_e)];
        end
        for ind_e=1:N2.num_elex_levels
        e_i_=[e_i_, N2.ev_i{ind_e}+N2.ev_0(ind_e)+N2.e_E(ind_e)];
        end

    if (i_ini<=6)
        En0=n0*e_i*n/n1 + n0*e_i_*n_/n1 + n0*k*T0 + 1.5*n0*k*T0 +...
        n0*NO.form_e*f + n0*N2.form_e*(1-f);
        Ep0=(En0+n0*k*T0)/rho0+0.5*v0^2;  
        % (E0+p0)/rho0+v0^2/2
    else
        En0=n0*e_i*n/n1 + n0*e_i_*n_/n1 + (f + (1-f)/2)*n0*k*T0 + ...
        1.5*n0*k*T0 + n0*NO.form_e*f + n0*N2.form_e*(1-f);
        Ep0=(En0+n0*k*T0)/rho0+0.5*v0^2;  
    end
    if i_rel==2
    disp('Conservation laws check behind ISW')
    check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 0);
    end
    
    %% REFL
    
    if (i_exch==2)||(i_exch==3)
    %taking into account chemical processes
    Reacs_keys={'Diss','Exch', 'VT', 'VV'};
    reacs_val={Diss, Exch, model_VT, model_VT};
    elseif i_exch==1
    %without 
    Reacs_keys={'Diss', 'VT', 'VV'};
    reacs_val={Diss,  model_VT, model_VT};
    else
        disp("Error, invalid value of i_exch");
    end

    kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
    if i_rel==2
    n0=sum(Y(end, 1:end-2),2);   % m-3
    v0=v0+v0_r-Y(end, end-1);   % m/s
    T0=Y(end, end);   % K
    rho0=n_N2(end)*N2.mass + n_NO(end)*NO.mass + n_O2(end)*O2.mass...
        + n_O(end)*O.mass + n_N(end)*N.mass + n_Ar(end)*Ar.mass;
    elseif i_rel==1
    n0=sum(Y(1:end-2));
    v0=v0+v0_r-Y(end-1);
    T0=Y(end);
        if (i_ini<=6)
            rho0=n0*((1-f)*N2.mass + f*NO.mass);
        else
            rho0=n0*((1-f)/2*N2.mass + (1-f)/2*Ar.mass + f*NO.mass);
        end
    end
    if (i_ini<=6)
    [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0, 1);
    else 
    [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0, f+(1-f)/2);
    end
    kinetics.n0=n0;
    kinetics.v0=v0;
    kinetics.T0=T0;
    kinetics.Delta=Delta;
    timewave=1000*1e-6;
    x_w=v0_r*timewave;
    xspan=[0 x_w]./Delta;
    y0_1=zeros(kinetics.num_eq+2, 1);
    if i_rel==2
    y0_1(1:end)=Y(end, :).*((1/n0)*n1);
    elseif i_rel==1
    n=density_f_exc(T0buf, n1*f, NO);
    y0_1(kinetics.index{IndexOfMolecules("NO")})=n;
    if (i_ini<=6)
        n_=density_f_exc(T0, n1*(1-f), N2);
        y0_1(kinetics.index{IndexOfMolecules("Ar")})=0;
    else
        n_=density_f_exc(T0, n1*(1-f)/2, N2);
        y0_1(kinetics.index{IndexOfMolecules("Ar")})=n1*(1-f)/2;
    end
    y0_1(kinetics.index{IndexOfMolecules("N2")})=n_;
    end
    y0_1(end-1)=v1;
    y0_1(end)=T1;
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
    if i_exch~=3  %if i_exch==3 then TvNO has no meaning
    Tv_NO_1 = NO.ev_i{1}(2)./(k*log(Y_1(:, ...
      kinetics.index{IndexOfMolecules("NO")}(1))./...
      Y_1(:,kinetics.index{IndexOfMolecules("NO")}(2))));
    end
    Tv_O2=O2.ev_i{1}(2)./(k*log(Y_1(:, ...
        kinetics.index{IndexOfMolecules("O2")}(1))./...
     Y_1(:,kinetics.index{IndexOfMolecules("O2")}(2))));
    Tv_N2=N2.ev_i{1}(2)./(k*log(Y_1(:, ...
        kinetics.index{IndexOfMolecules("N2")}(1))./...
     Y_1(:,kinetics.index{IndexOfMolecules("N2")}(2))));
    time_ms_1=X_1./v0_r*1e6;
    if i_rel==2
    rho0=n_NO(end)*NO.mass + n_O(end)*O.mass + n_N(end)*N.mass + ...
        n_O2(end)*O2.mass + n_N2(end)*N2.mass + n_Ar(end)*Ar.mass;
    end

    rhov0_1=rho0 * v0;                    % rho0*v0
    rhov2p0_1=rho0* v0^2 + n0*k*T0;      % rho0*v0^2+p0
    En0_1=n0*e_i*y0_1(kinetics.index{IndexOfMolecules("NO")})/n1 +...
    n0*e_i_*y0_1(kinetics.index{IndexOfMolecules("N2")})/n1 ...
    + k*T0*n_NO(end) + k*T0*n_N2(end)  + 1.5*n0*k*T0 +...
    n_NO(end)*NO.form_e + n_N(end)*N.form_e+n_O(end)*O.form_e +...
    n_N2(end) *N2.form_e;
    Ep0_1=(En0_1+n0*k*T0)/rho0+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
    disp('Conservation laws check behind RSW')
    check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 0);

    %This is where the output data is stored. 
    % They contain the evolution of temperatures, number densities,
    % and pressure between the SWs and behind the reflected SW
    if i_rel==2
    resSt.time=time_ms;
    resSt.T=T;
    if i_exch~=3
    resSt.Tv=Tv;
    end
    resSt.nO=n_O;
    resSt.nN=n_N;
    resSt.nNO=n_NO;
    resSt.nAr=n_Ar;
    resSt.nO2=n_O2;
    resSt.nN2=n_N2;
    resSt.p=p;
    dat(i_vibr,i_U,i_ini, i_exch)=resSt;
    end
    resSt_1.time=time_ms_1;
    resSt_1.T=T_1;
    if i_exch~=3
    resSt_1.TvNO=Tv_NO_1;
    end
    resSt_1.TvO2=Tv_O2;
    resSt_1.TvN2=Tv_N2;
    resSt_1.ni_NO=Y_1(:, kinetics.index{IndexOfMolecules("NO")});
    resSt_1.ni_O2=Y_1(:, kinetics.index{IndexOfMolecules("O2")});
    resSt_1.ni_N2=Y_1(:, kinetics.index{IndexOfMolecules("N2")});
    resSt_1.p=p_1;
    resSt_1.nNO=n_NO_1/Na;
    resSt_1.nO2=n_O2_1/Na;
    resSt_1.nN2=n_N2_1/Na;
    resSt_1.nO=n_O_1/Na;
    resSt_1.nN=n_O_1/Na;
    resSt_1.nAr=n_Ar_1/Na;
    dat1(i_vibr,i_U,i_ini,i_rel, i_exch)=resSt_1;
end
end
end
end
end

%%
%if you want to save your data in .mat file, uncomment following raws
%save(['..\data\NO_N2 Streicher experiment\NO_N2_betweenSWs_output.mat'], 'dat');
%save(['..\data\NO_N2 Streicher experiment\NO_N2_behindRSW_output.mat'], 'dat1');
rmpath('../src/')
rmpath('../data/')
toc       