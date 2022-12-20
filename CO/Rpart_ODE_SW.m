function out=Rpart_ODE_SW(x, y, n0, T0, v0, Delta, Prcl, Coll, ...
                        ind_Arr, ind_U, i_dis, ind_exc, setup) %#ok<INUSL>
% Right part function for ODE systems for shock waves (SW).
% x is the x coordinate; y is the vector of macroparameters;
% n0 is the characteristic density; T0 is the characteristic temperature;
% v0 is the characteristic velocity, Delta is the characteristic length;
% Prcl is the particles container; Coll is the collisions container;
% ind_Arr is the Arrhenius law constants indicator; 
% ind_U is the indicator of the non-equilibrium parameter U in 
% dissociation models; i_dis is the dissociation and kinetic scheme 
% indicator; ind_exc is the indicator of the electronic excitation 
% presence; setup is the different parameters for kinetic scheme container
% 19.12.2022 Maksim Melnik

ind_Aliat=0;
if i_dis>2
    ind_Aliat=1;
end
num_vibr_levels=Prcl.CO.num_vibr_levels(1);
num_COa=num_vibr_levels+1;
num_COA=num_COa+Prcl.CO.num_vibr_levels(2);
num_C=num_COA+Prcl.CO.num_vibr_levels(3);
num_O=num_C+1;
num_C2=num_O+1;
num_Ar=num_C2+1;
k=1.380648528e-23;
ni_b = y(1:num_vibr_levels);
nCOa_b = y(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
nCOA_b = y(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
nac_b = y(num_C);
nao_b = y(num_O);
nC2_b = y(num_C2);
naAr_b = y(num_Ar);

% working on universal Rci
N_a=6.02214076e23;          % Avogadro constant
switch ind_Arr
    case 1
        Diss.Arrhenius='Park';
        keys={'CO', 'C', 'O', 'Ar', 'C2'};
        val_A={2.3e20/N_a*1e-6, 3.4e20/N_a*1e-6, 3.4e20/N_a*1e-6, ...
                                    2.3e19/N_a*1e-6, 2.3e20/N_a*1e-6};
        val_n={-1, -1, -1, -1, -1};
        diss_Arrhenius_A_CO=containers.Map(keys, val_A);
        diss_Arrhenius_n_CO=containers.Map(keys, val_n);
        val_A={3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6, ...
                                    3.7e14/N_a*1e-6, 3.7e14/N_a*1e-6};
        val_n={0, 0, 0, 0, 0};
        diss_Arrhenius_A_C2=containers.Map(keys, val_A);
        diss_Arrhenius_n_C2=containers.Map(keys, val_n);
    case 2
        Diss.Arrhenius='Ibraguimova';
    case 3
        Diss.Arrhenius='McKenzie';
    case 4
        Diss.Arrhenius='Fairbairn';
end
mix.num=0;
CO=Prcl.CO;
CO.diss_Arrhenius_A=diss_Arrhenius_A_CO;
CO.diss_Arrhenius_n=diss_Arrhenius_n_CO;
Ps={mix, CO, Prcl.C, Prcl.O};
if setup.C2
    C2=Prcl.C2;
    C2.diss_Arrhenius_A=diss_Arrhenius_A_C2;
    C2.diss_Arrhenius_n=diss_Arrhenius_n_C2;
    Ps{length(Ps)+1}=C2;
end
if setup.f<1
    Ps{length(Ps)+1}=Prcl.Ar;
end
num=0;
index{1}=0;
for ind=2:length(Ps)
    num_states=sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    num=num+num_states;
    first=index{ind-1}(end)+1;
    index{ind}=first:first+num_states-1;
end
Ps{1}.num=num;

    % call of relaxation terms Rci
setupV2=setup;
setupV2.T0=T0;
setupV2.n0=n0;
setupV2.ind_Arr=ind_Arr;
setupV2.ind_U=ind_U;
setupV2.ind_Aliat=ind_Aliat;
Reacs_keys={'Diss', 'VT'};  % VE?
if ind_exc
    Reacs_keys=[Reacs_keys, 'VE'];
end
Diss.rec=true;              % Is recombination included?
switch i_dis
    case 1
        Diss.NEmodel='MT';
    case 2
        Diss.NEmodel='MT';
    case 3
        Diss.NEmodel='Aliat';
    case 4
        Diss.NEmodel='Savelev21';
end
switch ind_U
    case 2
        Diss.U='D/6k';
    case 3
        Diss.U='3T';
    case 4
        Diss.U='inf';
end
kinetics.Ps=Ps(2:end);
kinetics.num_Ps=length(kinetics.Ps);
kinetics.num_eq=num;
reacs_val={Diss, setup.model_VT};
if ind_exc
    reacs_val=[reacs_val, 'VE'];
end
kinetics.reactions=containers.Map(Reacs_keys, reacs_val);
kinetics.index=index(2:end);
kinetics.T0=T0;
kinetics.n0=n0;

y2=ni_b;
if ind_exc
    y2=[y2; nCOa_b; nCOA_b];
end
y2=[y2; nac_b; nao_b];
if setup.C2
    y2=[y2; nC2_b];
end
if setup.f<1
    y2=[y2; naAr_b];
end
y2=[y2; y(end-1:end)];
R2=Rci(y2, kinetics);

v_DN=y2(end-1);     % dimentionless gas velocity
T_DN=y2(end);
    % number densities equations
M2 = diag([ones(1, kinetics.num_eq)*v_DN 0 0]);
M2(1:end-2, end-1) = y2(1:end-2);
    % momentum equation
M2(end-1, 1:kinetics.num_eq) = T_DN;
rho=0;
for ind=1:kinetics.num_Ps
    rho = rho + kinetics.Ps{ind}.mass*sum(y2(kinetics.index{ind}));
end
M2(end-1, end-1) = rho*v_DN*v0^2/k/T0;
M2(end-1, end)   = sum(y2(1:end-2));
    % energy equation
nm=0;   % number density of molecules
for ind=1:kinetics.num_Ps
 if kinetics.Ps{ind}.fr_deg_c>3
  e_i=[];
  for ind_e=1:kinetics.Ps{ind}.num_elex_levels
   e_i=[e_i, kinetics.Ps{ind}.ev_i{ind_e}+kinetics.Ps{ind}.ev_0(ind_e)+...
       kinetics.Ps{ind}.e_E(ind_e)];
  end
  M2(end, kinetics.index{ind}) = (e_i +kinetics.Ps{ind}.form_e)/k/T0 ...
                                                                +2.5*T_DN;
  M2(end, end-1)=M2(end, end-1)+((e_i +kinetics.Ps{ind}.form_e)/k/T0)*...
      y2(kinetics.index{ind})+3.5*sum(y2(kinetics.index{ind}))*T_DN;
  nm=nm+sum(y2(kinetics.index{ind}));
 else
  M2(end, kinetics.index{ind}) = (kinetics.Ps{ind}.form_e ...
      + kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels))/k/T0 ...
                                                                +1.5*T_DN;
  M2(end, end-1)=M2(end, end-1)+2.5*sum(y2(kinetics.index{ind}))*T_DN+...
      (kinetics.Ps{ind}.e_E+kinetics.Ps{ind}.form_e)/k/T0*...
      y2(kinetics.index{ind});
 end
end
M2(end, end-1)=M2(end, end-1)/v_DN;
M2(end, end)=1.5*sum(y2(1:end-2))+nm;
M2sp=sparse(M2);
R2=[R2; 0; 0]*n0*Delta/v0;
out_temp=M2sp^(-1)*R2;
if ind_exc==0
    out_temp=[out_temp(1:68); zeros(120-68, 1); out_temp(69:end)];
end
out=out_temp;

end