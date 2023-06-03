function R2=Rci(y, kinetics)
% Universal function for relaxation terms R_{c\alpha i}.
% y is the vector of gas macroparameters on the current solving step;
% kinetics is the big structure with all kinetics.
% 09.12.2022 by Maksim Melnik.

T=y(end)*kinetics.T0;
n0=kinetics.n0;
R_VT_data2=zeros(kinetics.num_eq, 1);
R_VV_data=zeros(kinetics.num_eq, 1);
R_diss_data2=zeros(kinetics.num_eq, 1);
R_VE_data2=zeros(kinetics.num_eq, 1);
R_exch_data2=zeros(kinetics.num_eq, 1);
y2=y;


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
    for indM2=indM1:kinetics.num_Ps
     M2=kinetics.Ps{indM2};
     if M2.num_vibr_levels(1)>1
      i2=kinetics.index{indM2};
      for ind_e1=1:M1.num_elex_levels
       i1_e=i1(1+sum(M1.num_vibr_levels(1:ind_e1-1)):...
                                       sum(M1.num_vibr_levels(1:ind_e1)));
       for ind_e2=1:M2.num_elex_levels
        i2_e=i2(1+sum(M2.num_vibr_levels(1:ind_e2-1)):...
                                       sum(M2.num_vibr_levels(1:ind_e2))); 
        R_VV_data_temp=R_VV_old(M1, y2(i1_e), M2, y2(i2_e), ...
                            T, ind_e1, ind_e2, kinetics.reactions('VV'));
        R_VV_data(i1_e)=R_VV_data(i1_e)+sum(R_VV_data_temp, 2);
        if M2.name~=M1.name
         R_VV_data(i2_e)=R_VV_data(i2_e)+sum(R_VV_data_temp, 1)';
        end
       end
      end
     end
    end
   end
   
  end
  
  if M1.num_elex_levels>1
   if isKey(kinetics.reactions, 'VE')
    if M1.name=="CO"
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
    else
        str_w=strcat('VE is allowed only for CO, but ', M1.name, ...
                                        ' is electronicaly excited too.');
        warning(str_w)
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
   for indM2=1:kinetics.num_Ps
    M2=kinetics.Ps{indM2};
    i2=kinetics.index{indM2};
    coll2.ArrA=M1.diss_Arrhenius_A(M2.name);
    coll2.ArrN=M1.diss_Arrhenius_n(M2.name);
    y_diss=y2(i1);
    switch kinetics.reactions('Diss').NEmodel
     case {'MT', 'Aliat'}
      R_diss_data_temp=R_diss_Aliat_onoff(M1, y_diss, nP1, nP2, ...
           sum(y2(i2)), coll2, T, n0, kinetics.reactions('Diss').U, ...
                                    kinetics.reactions('Diss').NEmodel);
     case 'Savelev21'
       error("Savelev's diss model is still not implemented")
    end
    R_diss_data2(i1)=R_diss_data2(i1)+R_diss_data_temp;
   end
   R_diss_data2(iP1)=R_diss_data2(iP1)-sum(R_diss_data2(i1));
   R_diss_data2(iP2)=R_diss_data2(iP2)-sum(R_diss_data2(i1));
  end

   if isKey(kinetics.reactions, 'Exch') %exchange reactions
   if (M1.name=="NO") %first reaction NO+O->O2 + N
       %M1=NO
    indO=kinetics.index{2};
    indN=kinetics.index{3};
    indO2=kinetics.index{4};
    indNO=kinetics.index{1};
    N_a=6.02214076e23; 
    coll.ArrA(1)=8.4e12/N_a*1e-6; coll.ArrN(1)=0; %Arrhenius parameters from Park
    R_exch_temp=R_exch_NO_O__O2_N(M1, kinetics.Ps{4} , y2(indNO),...
        y2(indO), y2(indO2),  y2(indN), coll, T);
    %если я правильно понимаю для тех кто слева надо +, а для тех кто
    %справа -
    R_exch_data2(indNO)=R_exch_data2(indNO) + sum(R_exch_temp,1)'; 
    R_exch_data2(indO2)=R_exch_data2(indO2) - sum(R_exch_temp,2);
    R_exch_data2(indO)=R_exch_data2(indO) + sum(R_exch_data2(indNO));
    R_exch_data2(indN)=R_exch_data2(indN) - sum(R_exch_data2(indNO));
   end
   if (M1.name=="N2") %second reaction N2 + O -> NO + N
    indO=kinetics.index{2};
    indN=kinetics.index{3};
    indN2=kinetics.index{5};
    indNO=kinetics.index{1};
    N_a=6.02214076e23; 
    coll.ArrA(1)=6.4e17/N_a*1e-6; coll.ArrN(1)=-1; %Arrhenius parameters from Park
    R_exch_temp=R_exch_N2_O__NO_N(M1, kinetics.Ps{1} , y2(indN2),...
        y2(indO), y2(indNO),  y2(indN), coll, T);
    %если я правильно понимаю для тех кто слева надо +, а для тех кто
    %справа -
    R_exch_data2(indN2)=R_exch_data2(indN2) + sum(R_exch_temp,1)';
    R_exch_data2(indNO)=R_exch_data2(indNO) - sum(R_exch_temp,2);  
    R_exch_data2(indO)=R_exch_data2(indO) + sum(R_exch_data2(indN2));
    R_exch_data2(indN)=R_exch_data2(indN) - sum(R_exch_data2(indN2));
   end
 end
  
 end
end

R2=R_VT_data2+R_VV_data+R_diss_data2+R_VE_data2+R_exch_data2;

end