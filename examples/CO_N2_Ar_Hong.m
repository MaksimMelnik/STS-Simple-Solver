% The main function for the macroparameters calculation behind reflected SW
% for Streicher's experiment conditions in CO-Ar mixture.
%[1] J. Streicher, A. Krish, R. Hanson, High-temperature vibrational
%relaxation and decomposition of shock-heated nitric oxide: II. Nitrogen
%dilution from 1900 to 8200 K, Physics of Fluids 34 (11) (2022) 116123.
%doi:10.1063/5.0122787
%06.04.2023 Denis Kravchenko

tic
clearvars;
% constants
k=1.380649e-23;  % Boltzmann constant, J/K
Torr=133.322368;
Na=6.02214076e23;

%initialization of structures dat and dat1
tmp.time=0; tmp.T=0; tmp.Tv=0; tmp.nCO=0; tmp.nAr=0;
tmp.nN2=0; tmp.p=0; tmp.ni_CO = 0;
dat(2,4,12,3)=tmp;

tmp1.time=0; tmp1.T=0; tmp1.TvCO=0; tmp1.TvN2=0; tmp1.ni_CO=0;
tmp1.ni_N2=0; tmp1.p=0; tmp1.nCO=0; tmp1.nN2=0;
tmp1.nAr=0;
dat1(2,4,12,2,3)=tmp1;

clear tmp tmp1;
addpath('../src/')
addpath('../data/')
load('particles.mat', "CO", "N2", "Ar");
load('../data/reactions.mat'); %load reaction data

CO.num_elex_levels=1;       % CO electronic excitation
N2.num_elex_levels=1;
Ar.num_elex_levels=1;

% initial conditions
init_c= [ % f;  p0, Pa;   v0, m/s;   T0, K;   v0_1
    0.01 0.1 0.89 2190 1120 298 562       % 1% CO; 10% N2; 89% Ar
    0.01 0.05 0.94 2430 1101 298 573     % 1% CO; 5% N2; 94% Ar
    0.01 0.01 0.98 3410 1015 298 551    % 1% CO; 1% N2; 98% Ar
    ];

dbg = 1;

for i_ini=3 % [1 2 3]
    %choosing testcase

    for i_vibr=2 %[1 2]
        % choosing desired vibrational energy exchange model 1 for SSH; 2 for FHO

        for i_rel=2 %[1 2]
            % 1 -relaxation off; 2 - relaxation on

            fCO=init_c(i_ini, 1); %molar fraction of CO
            fN2=init_c(i_ini, 2); %molar fraction of N2
            fAr=init_c(i_ini, 3); %molar fraction of Ar
            fMol = fCO + fN2;
            assert(fCO+fN2+fAr == 1)
            p0=init_c(i_ini, 4); %initial pressure in shock tube
            v0=init_c(i_ini, 5);   % velocity of incident SW, m/s
            v0_i=v0;
            v0_r=init_c(i_ini, 7); %velocity of reflected SW, m/s
            T0=init_c(i_ini, 6); % initial temperature in shock tube, K
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

            %without exchange and diss-rec reactions, since this can be neglected
            Reacs_keys={'VT', 'VV'};
            reacs_val={model_VT, model_VT};
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
            n_boltz_CO=density_f_exc(T0, n1*fCO, CO);
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
            elseif i_rel==1
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
            % rho0*v0
            rhov2p0=rho0* v0^2 + n0*k*T0;
            % rho0*v0^2+p0
            e_i_CO=[];
            e_i_N2=[];
            for ind_e=1:CO.num_elex_levels
                e_i_CO=[e_i_CO, CO.ev_i{ind_e}+CO.ev_0(ind_e)+CO.e_E(ind_e)];
            end
            for ind_e=1:N2.num_elex_levels
                e_i_N2=[e_i_N2, N2.ev_i{ind_e}+N2.ev_0(ind_e)+N2.e_E(ind_e)];
            end
            En0=n0*e_i_CO*n_boltz_CO/n1 + n0*e_i_N2*n_boltz_N2/n1 + n0*k*T0 + 1.5*n0*k*T0;
            Ep0=(En0+n0*k*T0)/rho0+0.5*v0^2;
            % (E0+p0)/rho0+v0^2/2
        end
        if i_rel==2
            disp('Conservation laws check behind ISW')
            check_CL_SW([rhov0 rhov2p0 Ep0], Y, kinetics, 1);
        end

        %% REFL
        Reacs_keys={'VT', 'VV'};
        reacs_val={model_VT, model_VT};
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
        end
        [n1, v1, T1]=in_con_SW(n0, v0, T0, rho0, fMol);

        kinetics.n0=n0;
        kinetics.v0=v0;
        kinetics.T0=T0;
        kinetics.Delta=Delta;
        %time interval of calculation behind RSW
        timewave=400*1e-6;
        x_w=v0_r*timewave;
        xspan=[0 x_w]./Delta;
        y0_1=zeros(kinetics.num_eq+2, 1);

        %the vector of initial values тАЛтАЛfor modeling a reflected SW,
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
            n_boltz_N2=density_f_exc(T0, n1*fN2, N2);
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
        En0_1=n0*e_i_CO*y0_1(kinetics.index{IndexOfMolecules("CO")})/n1 +...
            n0*e_i_N2*y0_1(kinetics.index{IndexOfMolecules("N2")})/n1 ...
            + k*T0*n_CO(end) + k*T0*n_N2(end)  + 1.5*n0*k*T0;
        Ep0_1=(En0_1+n0*k*T0)/rho0+0.5*v0^2;       % (E0+p0)/rho0+v0^2/2
        disp('Conservation laws check behind RSW')
        check_CL_SW([rhov0_1 rhov2p0_1 Ep0_1], Y_1, kinetics, 1);

        %This is where the output data is stored.
        % They contain the evolution of temperatures, number densities,
        % and pressure between the SWs and behind the reflected SW

        if i_rel==2
            resSt.time=time_ms;
            resSt.T=T;
            resSt.Tv=ones(length(time_ms),1)*NaN;
            resSt.nCO=n_CO;
            resSt.nAr=n_Ar;
            resSt.nN2=n_N2;
            resSt.p=p;
            resSt.ni_CO = Y(:, kinetics.index{IndexOfMolecules("CO")});
            resSt_1.ni_N2=Y_1(:, kinetics.index{IndexOfMolecules("N2")});
            dat(i_vibr,i_ini)=resSt;
        end
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
        dat1(i_vibr,i_ini,i_rel)=resSt_1;
    end
end
figure
plot(time_ms_1/1e3, Tv_N2, time_ms_1/1e3, Tv_N2, 'linewidth', 1.5)

%%
%if you want to save your data in .mat file, uncomment following raws
% save('..\data\CO_N2 Streicher experiment\CO_N2_betweenSWs_output.mat', 'dat');
% save('..\data\CO_N2 Streicher experiment\CO_N2_behindRSW_output.mat', 'dat1');
rmpath('../src/')
rmpath('../data/')
toc