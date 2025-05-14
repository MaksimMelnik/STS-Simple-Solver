%The main function for the macroparameters calculation behind reflected SW
%for Hong's experiment conditions in CO-NO-Ar mixture.
%[1]He D. et al. Vibrational energy relaxation in shock-heated CO/N2/Ar mixtures
%The Journal of Chemical Physics. – 2024. – Т. 160. – №. 22.
%https://pubs.aip.org/aip/jcp/article-abstract/160/22/224308/3298204/Vibrational-energy-relaxation-in-shock-heated-CO
%20.04.2025 Semyon Zernov. This program inspired by Denis Kravchenko's code



tic
clearvars;
% constants
k=1.380649e-23;  % Boltzmann constant, J/K
Torr=133.322368;
Na=6.02214076e23;

dbg = 1;
vt_only = false;

tmp.time=0;
tmp.T=0;
tmp.p=0;
tmp.TvCO=0; tmp.TvN2=0;
tmp.ni_CO=0; tmp.ni_N2=0;
tmp.nCO=0; tmp.nN2=0; tmp.nAr=0;
result(3, 14, 3, 2)=tmp;

clear tmp;
addpath('../src/')
addpath('../data/')
load('particles.mat', "CO", "N2", "Ar");
load('../data/reactions.mat'); %load reaction data

CO.num_elex_levels=1;       % CO electronic excitation
N2.num_elex_levels=1;
Ar.num_elex_levels=1;

times_table = [
    1000 700 500 400 400 300 700 600 1400 1000 NaN NaN NaN NaN
    1000 1200 600 400 500 400 1000 1000 700 500 NaN NaN NaN NaN
    1000 1200 700 500 500 400 300 300 1500 500 700 500 400 300
    ];

init_c_1 = [
    3180 1030.23 297.35 526.83
    2630 1074.14 297.25 544
    2210 1113.68 297.75 559.7
    1820 1161.05 298.15 579.65
    1420 1237.85 298.35 609.65
    1620 1205.56 298.35 596.58
    2960 1057.55 297.55 537.54
    2190 1120.47 297.75 562.4
    4850 931.03 297.35 488.7
    3440 1012.03 297.35 519.76
    ];
Tvib_exp1 = [
    705 780 615 814 742 860 525 785 550 700
    ];

init_c_2 = [
    3770 988.27 297.15 525.68
    5140 911.61 297.15 494.47
    2690 1073.61 297.45 561.21
    2240 1116.55 297.45 579.31
    1850 1158.64 296.85 597.1
    1620 1201.31 297.15 615.41
    4530 936.7 297.55 504.69
    4130 967.9 297.65 517.42
    2910 1050.04 297.95 551.43
    2430 1100.99 297.95 572.81
    ];
Tvib_exp2 = [
    620 450 700 800 705 950 480 648 746 860
    ];

init_c_3 = [
    5040 911.67 297.55 506.16
    5610 884.49 298.15 494.76
    4010 968.39 298.35 530.64
    3410 1014.96 298.45 550.91
    3020 1045.65 298.55 564.39
    2650 1079.84 298.55 579.49
    2210 1124.11 297.75 599.07
    1910 1162.02 298.15 616.12
    6230 862.02 297.95 485.29
    3770 991.28 297.65 540.44
    4530 939.33 297.35 517.94
    3230 1026.7 297.75  555.93
    2440 1094.89  298.25 586.13
    2050 1145.84 298.35 608.89
    ];
Tvib_exp3 = [
    525 600 500 735 730 750 650 905 460 730 640 690 640 890
    ];

mixture = [
    % fCO fN2 fAr
    0.01 0.1 0.89  % 1% CO; 10% N2; 89% Ar
    0.01 0.05 0.94 % 1% CO; 5% N2; 94% Ar
    0.01 0.01 0.98 % 1% CO; 1% N2; 98% Ar
    ];

for i_mixture=1:3
    fprintf('\n\t\t\t\t\tСмесь %d\n', i_mixture);
    fCO=mixture(i_mixture, 1); % molar fraction of CO
    fN2=mixture(i_mixture, 2); % molar fraction of N2
    fAr=mixture(i_mixture, 3); % molar fraction of Ar
    fMol = fCO + fN2;
    assert(fCO+fN2+fAr == 1)
    switch i_mixture
        case 1
            init_c = init_c_1;
            Tvib_exp = Tvib_exp1;
        case 2
            init_c = init_c_2;
            Tvib_exp = Tvib_exp2;
        case 3
            init_c = init_c_3;
            Tvib_exp = Tvib_exp3;
    end

    for i_ini=1:length(init_c)
        fprintf('\t\t\tЭксперимент %d\n', i_ini);

        for i_vibr=1:2
            % choosing desired vibrational energy exchange model 1 for SSH; 2 for FHO
            switch i_vibr
                case 1
                    kinetics_model = 'SSH';
                case 2
                    kinetics_model = 'FHO';
            end
            fprintf('\t\t%s\n', kinetics_model);
            for i_rel=1:3
                switch i_rel
                    case 1
                        model_name = 'Модель замороженной релаксации';
                    case 2
                        model_name = 'Модель частичной релаксации';
                    case 3
                        model_name = 'Верификационный метод';
                end
                fprintf('\t%s\n', model_name);
                % 1 - relaxation off;
                % 2 - relaxation on;
                % 3 - start with Tvib at t = 0


                p0=init_c(i_ini, 1); %initial pressure in shock tube
                v0=init_c(i_ini, 2);   % velocity of incident SW, m/s
                v0_i=v0;
                v0_r=init_c(i_ini, 4); %velocity of reflected SW, m/s
                T0=init_c(i_ini, 3); % initial temperature in shock tube, K
                T0buf=T0; %buffer variable for initial temperature
                n0=p0/(k*T0);   % initial number density, m-3
                n0buf=n0; %buffer variable for initial number density
                rho0=n0*(fN2*N2.mass + fAr*Ar.mass + fCO*CO.mass);

                [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0, fMol);
                %dimensionless variables

                sigma0 = pi*CO.diameter^2;
                Delta = 1 / sqrt(2) / n0 / sigma0; %free path length
                Ps = {CO, N2, Ar};

                switch i_vibr
                    case 1
                        model_VT='SSH';
                    case 2
                        model_VT='FHO';
                end

                if vt_only
                    Reacs_keys={'VT'};
                    reacs_val={model_VT};
                else
                    % without exchange and diss-rec reactions, since this can be neglected
                    Reacs_keys={'VT', 'VV'};
                    reacs_val={model_VT, model_VT};
                end

                kinetics.Ps = Ps;
                kinetics.num_Ps=length(kinetics.Ps);
                kinetics.index = indexes_for_Ps(Ps);
                kinetics.num_eq = kinetics.index{end}(end);
                kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
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
                deltas=20e-3; %distance from the transducer to the endwall
                timewave=(v0 + v0_r)/(v0*v1)*deltas/v0_r; %time between passing the SWs
                x_w=v0*timewave;
                xspan=[0 x_w]./Delta;

                %vector of initial values
                y0=zeros(kinetics.num_eq+2, 1);
                y0(end-1)=v1;
                y0(end)=T1;
                %initial distribution of vibrational level populations - Boltzmann
                %distribution for CO molecules
                n_boltz_CO=density_f_exc(T0, n1*fCO, CO);  % <---- Tvib
                y0(kinetics.index{IndexOfMolecules("CO")})=n_boltz_CO;

                %and for N2 molecules in accordance with their mole fractions
                n_boltz_N2=density_f_exc(T0, n1*fN2, N2);
                y0(kinetics.index{IndexOfMolecules("N2")})=n_boltz_N2;
                y0(kinetics.index{IndexOfMolecules("Ar")})=n1*fAr;
                if dbg == 0
                    % great for an accurate simulation
                    options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-18, ...
                        'NonNegative', 1:kinetics.num_eq+2);
                else
                    % enough for debugging
                    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                        'NonNegative', 1:kinetics.num_eq+2);
                end
                if i_rel==2 %if relaxation between SWs on
                    [X, Y]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics), xspan, ...
                        y0, options_s);
                    X=X*Delta;
                    Y(:, 1:end-2)=Y(:, 1:end-2)*n0;
                    Y(:, end-1)= Y(:, end-1)*v0 ;
                    Y(:, end)=Y(:, end)*T0;
                    n_N2 = sum(Y(:,kinetics.index{IndexOfMolecules("N2")}),2);
                    n_Ar=Y(:, kinetics.index{IndexOfMolecules("Ar")});
                    n_CO=sum(Y(:, kinetics.index{IndexOfMolecules("CO")}), 2);
                    T=Y(:, end);
                    p=(n_CO +n_Ar + n_N2).*k.*T;
                    time_ms=X./v0*1e6;
                elseif i_rel==1 || i_rel==3
                    %if relaxation between SWs is off, than using the R-H relation
                    Y=y0*n0;
                    Y(end)=y0(end)*T0;
                    T=Y(end);
                    Y(end-1)=v0*y0(end-1);
                    Y=Y';
                    n_Ar=Y(kinetics.index{IndexOfMolecules("Ar")});
                    n_CO=Y(kinetics.index{IndexOfMolecules("CO")});
                    n_N2=Y(kinetics.index{IndexOfMolecules("N2")});
                end

                rhov0=rho0 * v0;
                rhov2p0=rho0* v0^2 + n0*k*T0;

                e_i_CO=[];
                e_i_N2=[];
                for ind_e=1:CO.num_elex_levels
                    e_i_CO=[e_i_CO, CO.ev_i{ind_e}+CO.ev_0(ind_e)+CO.e_E(ind_e)];
                end
                for ind_e=1:N2.num_elex_levels
                    e_i_N2=[e_i_N2, N2.ev_i{ind_e}+N2.ev_0(ind_e)+N2.e_E(ind_e)];
                end
                En0=n0*e_i_CO*n_boltz_CO/n1 + n0*e_i_N2*n_boltz_N2/n1 + ...
                    fMol*n0*k*T0 + ...
                    1.5*n0*k*T0 + n0*CO.form_e*fCO + n0*N2.form_e*fN2;
                Ep0=(En0+n0*k*T0)/rho0+0.5*v0^2;
                % (E0+p0)/rho0+v0^2/2
                if i_rel==2
                    disp('Conservation laws check behind ISW')
                    check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 1);
                end

                %% REFL
                if vt_only
                    Reacs_keys={'VT'};
                    reacs_val={model_VT};
                else
                    % without exchange and diss-rec reactions, since this can be neglected
                    Reacs_keys={'VT', 'VV'};
                    reacs_val={model_VT, model_VT};
                end
                kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
                if i_rel==2
                    n0=sum(Y(end, 1:end-2),2);   % m-3
                    v0=v0+v0_r-Y(end, end-1);   % m/s
                    T0=Y(end, end);   % K
                    rho0=n_N2(end)*N2.mass + n_CO(end)*CO.mass + n_Ar(end)*Ar.mass;
                elseif i_rel==1
                    n0=sum(Y(1:end-2));
                    v0=v0+v0_r-Y(end-1);
                    T0=Y(end);
                    rho0=n0*(fN2*N2.mass + fAr*Ar.mass + fCO*CO.mass);
                elseif i_rel==3
                    n0=sum(Y(1:end-2));
                    v0=v0+v0_r-Y(end-1);
                    T0=Y(end);
                    rho0=n0*(fN2*N2.mass + fAr*Ar.mass + fCO*CO.mass);
                end
                [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0, fMol);

                kinetics.n0=n0;
                kinetics.v0=v0;
                kinetics.T0=T0;
                kinetics.Delta=Delta;

                %time interval of calculation behind RSW
                timewave=times_table(i_mixture,i_ini)*1e-6;
                x_w=v0_r*timewave;
                xspan=[0 x_w]./Delta;
                y0_1=zeros(kinetics.num_eq+2, 1);

                %the vector of initial values for modeling a reflected SW,
                %in the case of taking into account relaxation, is taken from the
                %vector obtained when solving the problem of an incident SW,
                %but taking into account new dimensionaless with respect
                %to the new number density
                if i_rel==2
                    y0_1(1:end)=Y(end, :).*((1/n0)*n1);
                    % in case without intermediate relaxation it's filled in the
                    % same way as y0 before an incident SW
                elseif i_rel==1
                    n_boltz_CO=density_f_exc(T0buf, n1*fCO, CO);
                    y0_1(kinetics.index{IndexOfMolecules("CO")})=n_boltz_CO;
                    y0_1(kinetics.index{IndexOfMolecules("Ar")})=n1*fAr;
                    n_boltz_N2=density_f_exc(T0buf, n1*fN2, N2);
                    y0_1(kinetics.index{IndexOfMolecules("N2")})=n_boltz_N2;
                elseif i_rel==3
                    n_boltz_CO=density_f_exc(Tvib_exp(i_ini), n1*fCO, CO);
                    y0_1(kinetics.index{IndexOfMolecules("CO")})=n_boltz_CO;
                    y0_1(kinetics.index{IndexOfMolecules("Ar")})=n1*fAr;
                    n_boltz_N2=density_f_exc(Tvib_exp(i_ini), n1*fN2, N2);
                    y0_1(kinetics.index{IndexOfMolecules("N2")})=n_boltz_N2;
                end
                y0_1(end-1)=v1;
                y0_1(end)=T1;
                if dbg == 0
                    % great for an accurate simulation
                    options_s = odeset('RelTol', 3e-14, 'AbsTol', 1e-25, ...
                        'NonNegative', 1:kinetics.num_eq+2);
                else
                    % enough for debugging
                    options_s = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, ...
                        'NonNegative', 1:kinetics.num_eq+2);
                end
                [X_1, Y_1]=ode15s(@(t, y) Rpart_ODE_SW(t, y, kinetics),...
                    xspan, y0_1, options_s);


                X_1=X_1*Delta;
                Y_1(:, 1:end-2)=Y_1(:, 1:end-2)*n0;
                Y_1(:, end-1)=Y_1(:, end-1)*v0;
                Y_1(:, end)=Y_1(:, end)*T0;
                T_1=Y_1(:, end);
                n_Ar_1=Y_1(:, kinetics.index{IndexOfMolecules("Ar")});
                n_CO_1=sum(Y_1(:, kinetics.index{IndexOfMolecules("CO")}), 2);
                n_N2_1=sum(Y_1(:, kinetics.index{IndexOfMolecules("N2")}), 2);
                p_1=(n_CO_1 + n_Ar_1 + n_N2_1)*k.*T_1 /Torr;
                Tv_CO=CO.ev_i{1}(2)./(k*log(Y_1(:, ...
                    kinetics.index{IndexOfMolecules("CO")}(1))./...
                    Y_1(:,kinetics.index{IndexOfMolecules("CO")}(2))));
                Tv_N2=N2.ev_i{1}(2)./(k*log(Y_1(:, ...
                    kinetics.index{IndexOfMolecules("N2")}(1))./...
                    Y_1(:,kinetics.index{IndexOfMolecules("N2")}(2))));
                time_ms_1=X_1./v0_r*1e6;
                if i_rel==2
                    rho0=n_CO(end)*CO.mass + n_N2(end)*N2.mass + n_Ar(end)*Ar.mass;
                end

                rhov0_1=rho0 * v0;                    % rho0*v0
                rhov2p0_1=rho0* v0^2 + n0*k*T0;      % rho0*v0^2+p0
                En0_1=n0*e_i_CO*y0_1(kinetics.index{IndexOfMolecules("CO")})/n1 + ...
                    n0*e_i_N2*y0_1(kinetics.index{IndexOfMolecules("N2")})/n1 ...
                    + k*T0*n_CO(end) + k*T0*n_N2(end) + 1.5*n0*k*T0 + ...
                    n_CO(end)*CO.form_e + n_N2(end)*N2.form_e;
                Ep0_1=(En0_1+n0*k*T0)/rho0+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
                disp('Conservation laws check behind RSW')
                check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 1);

                resSt_1.time=time_ms_1;
                resSt_1.T=T_1;
                resSt_1.TvCO=Tv_CO;
                resSt_1.TvN2=Tv_N2;
                resSt_1.ni_CO=Y_1(:, kinetics.index{IndexOfMolecules("CO")});
                resSt_1.ni_N2=Y_1(:, kinetics.index{IndexOfMolecules("N2")});
                resSt_1.p=p_1;
                resSt_1.nCO=n_CO_1/Na;
                resSt_1.nN2=n_N2_1/Na;
                resSt_1.nAr=n_Ar_1/Na;
                result(i_mixture,i_ini,i_rel,i_vibr)=resSt_1;
            end
        end
    end
end

if vt_only
    save('..\data\CO_N2_Ar Hong experiment\CO_N2_Ar_behindRSW_output_VT.mat', 'result');
else
    save('..\data\CO_N2_Ar Hong experiment\CO_N2_Ar_behindRSW_output_VT_VV.mat', 'result');
end
rmpath('../src/')
rmpath('../data/')
toc