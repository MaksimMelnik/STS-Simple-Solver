function R=Rci(y, Prcl, setup, ind_exc, kinetics)
% Universal function for relaxation terms R_{c\alpha i}.
% y is the vector of gas macroparameters on the current solving step;
% Prcl is the container for all particles; Coll is the container
% for all collision parameters; setup is the problem setup;
% ind_exc indicates if CO is excited; i_dis is the dissociation model 
% indicator; ind_Arr is the Arrenius constants indicator;
% ind_U is the U parameter indicator in dissociation models;
% ind_Aliat is the Aliat or MT switcher.
%  Test variables: Ps is the particles array, kinetics is the big 
% structure with all kinetics.
% 09.12.2022 by Maksim Melnik.

num_vibr_levels=Prcl.CO.num_vibr_levels(1);
num_COa=num_vibr_levels+1;
num_COA=num_COa+Prcl.CO.num_vibr_levels(2);
num_C=num_COA+Prcl.CO.num_vibr_levels(3);
num_O=num_C+1;
num_C2=num_O+1;
num_Ar=num_C2+1;
R=zeros(num_Ar, 1);
ni_b = y(1:num_vibr_levels);
nCOa_b = y(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
nCOA_b = y(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
nac_b = y(num_C);
nao_b = y(num_O);
nC2_b = y(num_C2);
naAr_b = y(num_Ar);

% the new format
T=y(end)*setup.T0;
n0=setup.n0;
R_VT_data2=zeros(kinetics.num_eq, 1);
R_diss_data2=zeros(kinetics.num_eq, 1);
R_VE_data2=zeros(kinetics.num_eq, 1);
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
%     disp('Ar is here')
end
y2=[y2; y(end-1:end)];
if isKey(kinetics.reactions, 'Exch')
    error('Exchange reactions are still not implemented.')
end
for indM1=1:kinetics.num_Ps     % considering each particle
 M1=kinetics.Ps{indM1};
 i1=kinetics.index{indM1};      % pointer on ni of M1
 if M1.fr_deg_c>3

  if M1.num_vibr_levels(1)>1     % should we consider vibr processes?
   if isKey(kinetics.reactions, 'VT')
    for indM2=1:kinetics.num_Ps
     M2=kinetics.Ps{indM2};
     i2=kinetics.index{indM2};
     for ind_e=1:M1.num_elex_levels
      i1_e=i1(1+sum(M1.num_vibr_levels(1:ind_e-1)):...
                                        sum(M1.num_vibr_levels(1:ind_e)));
      R_VT_data_temp=R_VT_old(M1, y2(i1_e), M2, ...
                        sum(y2(i2)), T, ind_e, kinetics.reactions('VT'));
      R_VT_data2(i1_e)=R_VT_data2(i1_e)+R_VT_data_temp;
     end
    end
   end
   
   if isKey(kinetics.reactions, 'VV')
       error('VV is still not implemented.')
   end
  end
  
  if M1.num_elex_levels>1
   if isKey(kinetics.reactions, 'VE')
%     error('VE is still not implemented.')
    i1_e1=i1(1:M1.num_vibr_levels(1));
    np=0;
    for indM2=1:kinetics.num_Ps
     np=np+sum(y2(kinetics.index{indM2}));
    end
    for ind_e=2:M1.num_elex_levels
     i1_e2=i1(1+sum(M1.num_vibr_levels(1:ind_e-1)):...
                                        sum(M1.num_vibr_levels(1:ind_e)));
     R_VE_temp=R_VE_m(M1, y2(i1_e1), y2(i1_e2), np, T, ind_e);
     R_VE_data2(i1_e1)=R_VE_data2(i1_e1)+sum(R_VE_temp, 2);
     R_VE_data2(i1_e2)=R_VE_data2(i1_e2)-sum(R_VE_temp, 1)';
    end
   end
  end
  
  if isKey(kinetics.reactions, 'Diss')
   for indM3=1:kinetics.num_Ps     % finding indexes of diss parts of M1
    if kinetics.Ps{indM3}.name==M1.diss_parts(1)
     indP1=indM3;
    end
    if kinetics.Ps{indM3}.name==M1.diss_parts(2)
     indP2=indM3;
    end
   end
   iP1=kinetics.index{indP1};
   iP2=kinetics.index{indP2};
   nP1=sum(y2(iP1));
   nP2=sum(y2(iP2));
%    что-то с диссоциацией С2 всё равно не так
%    disp(kinetics.Ps{indM1}.name)
   for indM2=1:kinetics.num_Ps
    M2=kinetics.Ps{indM2};
    i2=kinetics.index{indM2};
    coll2.ArrA=M1.diss_Arrhenius_A(M2.name);
    coll2.ArrN=M1.diss_Arrhenius_n(M2.name);
    y_diss=y2(i1);
    if kinetics.Ps{indM1}.name=="CO"
       if ind_exc==0
           y_diss=[y_diss; zeros(120-68, 1)];
       end
    end
    switch kinetics.reactions('Diss').NEmodel
     case {'MT', 'Aliat'}
      R_diss_data_temp=R_diss_Aliat_onoff(M1, y_diss, nP1, nP2, ...
           sum(y2(i2)), coll2, T, n0, setup.ind_Arr, setup.ind_U, ...
                                                        setup.ind_Aliat);
       if kinetics.Ps{indM1}.name=="CO"
        if ind_exc==0
           R_diss_data_temp=R_diss_data_temp(1:68);
        end
       end
     case 'Savelev21'
       error("Savelev's diss model is still not implemented")
    end
    R_diss_data2(i1)=R_diss_data2(i1)+R_diss_data_temp;
   end
   R_diss_data2(iP1)=R_diss_data2(iP1)-sum(R_diss_data2(i1));
   R_diss_data2(iP2)=R_diss_data2(iP2)-sum(R_diss_data2(i1));
  end
  
 end
end
R2=R_VT_data2+R_diss_data2+R_VE_data2;

R(1:num_vibr_levels)=R2(1:num_vibr_levels);
if ind_exc
    R(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1)=...
        R2(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
    R(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1)=...
        R2(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
end
R(num_C)=R2(kinetics.index{2});
R(num_O)=R2(kinetics.index{3});
if setup.C2
    R(num_C2)=R2(kinetics.index{4});
end
end