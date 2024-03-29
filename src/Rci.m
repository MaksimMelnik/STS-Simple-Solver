function [R, Qin] = Rci(y, kinetics)
% Universal function for relaxation terms R_{c\alpha i}.
% R is the relaxation term Rci, Qin is the energy flux.
% y is the vector of gas macroparameters on the current solving step;
% kinetics is the big structure with all kinetics.
% 09.12.2022 by Maksim Melnik.

T = y(end) * kinetics.T0;
n0 = kinetics.n0;
R_VT_data   = zeros(kinetics.num_eq, 1);
R_VV_data   = zeros(kinetics.num_eq, 1);
R_diss_data = zeros(kinetics.num_eq, 1);
R_VE_data   = zeros(kinetics.num_eq, 1);
R_exch_data = zeros(kinetics.num_eq, 1);
R_wall_data = zeros(kinetics.num_eq, 1);
R_ET_data   = zeros(kinetics.num_eq, 1);
Qin = 0;

for indM1 = 1:kinetics.num_Ps   % considering each particle
 M1 = kinetics.Ps{indM1};
 i1 = kinetics.index{indM1};    % pointer on ni of M1
 if M1.fr_deg_c > 3
      
  if isKey(kinetics.reactions, 'VT')
   for indM2 = 1:kinetics.num_Ps
    M2 = kinetics.Ps{indM2};
    i2 = kinetics.index{indM2};
    for ind_e = 1:M1.num_elex_levels
     if M1.num_vibr_levels(ind_e) > 1
      i1_e = i1(1+sum(M1.num_vibr_levels(1:ind_e-1)) : ...
                                        sum(M1.num_vibr_levels(1:ind_e)));
      [R_VT_data_temp, Q_VT] = R_VT(M1, y(i1_e), M2, ...
                        sum(y(i2)), T, ind_e, kinetics.reactions('VT'));
      R_VT_data(i1_e) = R_VT_data(i1_e) + R_VT_data_temp;
      Qin = Qin + Q_VT;
     end
    end
   end
  end
   
  if isKey(kinetics.reactions, 'VV')
   for ind_e1 = 1:M1.num_elex_levels
    if M1.num_vibr_levels(ind_e) > 1
     for indM2 = indM1:kinetics.num_Ps
      M2 = kinetics.Ps{indM2};
      for ind_e2 = 1:M2.num_elex_levels
       if M2.num_vibr_levels(ind_e2) > 1
        i2 = kinetics.index{indM2};
        i1_e = i1(1+sum(M1.num_vibr_levels(1:ind_e1-1)) : ...
                                       sum(M1.num_vibr_levels(1:ind_e1)));
        i2_e = i2(1+sum(M2.num_vibr_levels(1:ind_e2-1)) : ...
                                       sum(M2.num_vibr_levels(1:ind_e2))); 
        [R_VV_data_temp, Q_VV] = R_VV(M1, y(i1_e), M2, y(i2_e), ...
                             T, ind_e1, ind_e2, kinetics.reactions('VV'));
        R_VV_data(i1_e) = R_VV_data(i1_e) + sum(R_VV_data_temp, 2);
        Qin = Qin + Q_VV;
        if M2.name ~= M1.name
         R_VV_data(i2_e) = R_VV_data(i2_e) + sum(R_VV_data_temp, 1)';
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
      np=np+sum(y(kinetics.index{indM2}));
     end
     for ind_e=2:M1.num_elex_levels
      i1_e2=i1(1+sum(M1.num_vibr_levels(1:ind_e-1)):...
                                        sum(M1.num_vibr_levels(1:ind_e)));
      R_VE_temp=R_VE_m(M1, y(i1_e1), y(i1_e2), np, T, ind_e);
      R_VE_data(i1_e1)=R_VE_data(i1_e1)+sum(R_VE_temp, 2);
      R_VE_data(i1_e2)=R_VE_data(i1_e2)-sum(R_VE_temp, 1)';
     end
    else
        str_w=strcat('VE is allowed only for CO, but ', M1.name, ...
                                        ' is electronicaly excited too.');
        error(str_w)
    end
   end
  end
  
  if isKey(kinetics.reactions, 'Diss')
   for indM3 = 1:kinetics.num_Ps     % finding indexes of diss parts of M1
    if kinetics.Ps{indM3}.name == M1.diss_parts(1)
     indP1=indM3;
    end
    if kinetics.Ps{indM3}.name==M1.diss_parts(2)
     indP2=indM3;
    end
   end
   iP1=kinetics.index{indP1};
   iP2=kinetics.index{indP2};
   nP1=sum(y(iP1));
   nP2=sum(y(iP2));
   for indM2=1:kinetics.num_Ps
    M2=kinetics.Ps{indM2};
    i2=kinetics.index{indM2};
    coll2.ArrA=M1.diss_Arrhenius_A(M2.name);
    coll2.ArrN=M1.diss_Arrhenius_n(M2.name);
    y_diss=y(i1);
    [R_diss_data_temp, Q_diss] = R_diss(M1, y_diss, nP1, nP2, ...
          sum(y(i2)), coll2, T, n0, kinetics.reactions('Diss').U, ...
                                    kinetics.reactions('Diss').NEmodel);
    R_diss_data(i1) = R_diss_data(i1) + R_diss_data_temp;
    Qin = Qin + Q_diss;
   end
   R_diss_data(iP1) = R_diss_data(iP1) - sum(R_diss_data(i1));
   R_diss_data(iP2) = R_diss_data(iP2) - sum(R_diss_data(i1));
  end
  
  if isKey(kinetics.reactions, 'Wall')
   if M1.num_vibr_levels(1)>1
    if max(isKey(kinetics.reactions, {'VT', 'VV'}))
     [R_VT_wall_data_temp, Qwall] = R_VT_wall(M1, y(i1), T, kinetics);
     R_wall_data(i1) = R_wall_data(i1) + R_VT_wall_data_temp/kinetics.n0;
     Qin = Qin + Qwall/kinetics.n0;
    end
   end
   if isKey(kinetics.reactions, 'Rec_wall')
    if M1.sigma == 2 && M1.form_e == 0
     for indM3 = 1:kinetics.num_Ps   % finding indexes of diss parts of M1
      if kinetics.Ps{indM3}.name == M1.diss_parts(1)
       indP1 = indM3;
       break
      end
     end
     iP1 = kinetics.index{indP1};
     nP1 = y(iP1(1));
     [R_rec_wall_data_temp, Q_rec_wall] = ...
                        R_rec_wall(kinetics.Ps{indP1}, nP1, T, kinetics);
	 R_wall_data(iP1(1)) = R_wall_data(iP1(1)) + ...
                                        R_rec_wall_data_temp/kinetics.n0;
	 R_wall_data(i1(1)) = R_wall_data(i1(1)) - ...
                                   0.5 * R_rec_wall_data_temp/kinetics.n0;
	 Qin = Qin + Q_rec_wall/kinetics.n0;
    end
   end
   if isKey(kinetics.reactions, 'ET')
    if M1.num_elex_levels > 1
     [R_ET_data_temp, Q_ET] = R_ET_diff_wall(M1, y, T, kinetics);
     R_ET_data(i1) = R_ET_data(i1) + R_ET_data_temp/n0^2;
     Qin = Qin + Q_ET/n0^2;
    end
   end
  end
  
 end
end

if isKey(kinetics.reactions, 'Exch') % exchange reactions universal attempt
 IndexOfMolecules = kinetics.IndexOfMolecules;
 exch_reactions = kinetics.reactions("Exch");
 for ind_exch = 1:length(exch_reactions)
  reaction = exch_reactions(ind_exch);
  IOM_M1 = IndexOfMolecules(reaction.particles(1));
  IOM_M2 = IndexOfMolecules(reaction.particles(2));
  IOM_M3 = IndexOfMolecules(reaction.particles(3));
  IOM_M4 = IndexOfMolecules(reaction.particles(4));
  indM1 = kinetics.index{IOM_M1};
  indM2 = kinetics.index{IOM_M2};
  indM3 = kinetics.index{IOM_M3};
  indM4 = kinetics.index{IOM_M4};
  switch length(reaction.particles)
   case 4
%     indM1 = indM1(1:kinetics.Ps{IOM_M1}.num_vibr_levels(1));
    [R_exch_temp, Q_exch] = R_exch({kinetics.Ps{IOM_M1}, ...
       kinetics.Ps{IOM_M2}, kinetics.Ps{IOM_M3}, kinetics.Ps{IOM_M4}}, ...
                     y(indM1), y(indM2), y(indM3), y(indM4), T, reaction);
    R_exch_data(indM1) = R_exch_data(indM1) + sum(R_exch_temp, [2 3 4]);
    R_exch_data(indM2) = R_exch_data(indM2) + sum(R_exch_temp, [1 3 4])';
    R_exch_data(indM3) = R_exch_data(indM3) ...
                            - reshape(sum(R_exch_temp, [1 2 4]), [], 1);
    R_exch_data(indM4) = R_exch_data(indM4) ...
                            - reshape(sum(R_exch_temp, [1 2 3]), [], 1);
   case 5
    IOM_M5 = IndexOfMolecules(reaction.particles(5));
    indM5  = kinetics.index{IOM_M5};
    [R_exch_temp, Q_exch] = R_exch_23(...
       {kinetics.Ps{IOM_M1}, kinetics.Ps{IOM_M2}, kinetics.Ps{IOM_M3}, ...
        kinetics.Ps{IOM_M4}, kinetics.Ps{IOM_M5}}, ...
       y(indM1), y(indM2), y(indM3), y(indM4), y(indM5), T, reaction, n0);
    R_exch_data(indM1) = R_exch_data(indM1) + sum(R_exch_temp, [2, 3]);
    R_exch_data(indM2) = R_exch_data(indM2) + sum(R_exch_temp, [1, 3])';
    R_exch_data(indM3) = R_exch_data(indM3) ...
                            - reshape(sum(R_exch_temp, [1, 2]), [], 1);
    R_exch_data(indM4) = R_exch_data(indM4) ...
                            - reshape(sum(R_exch_temp, [1 2 3 5]), [], 1);
    R_exch_data(indM5) = R_exch_data(indM5) ...
                            - reshape(sum(R_exch_temp, [1 2 3 4]), [], 1);
   otherwise
       error(["Reactions with current number of particles are not " ...
                                                    "implemented yet"])
  end
  Qin = Qin + Q_exch;
 end
end

R = R_VT_data + R_VV_data + R_diss_data + R_VE_data + ... 
                                    R_exch_data + R_wall_data + R_ET_data;
end