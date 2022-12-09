function R=Rci(y, T0, n0, Prcl, Coll, setup, ind_exc, i_dis, ind_Arr, ...
                                                        ind_U, ind_Aliat)
% Universal function for relaxation terms R_{c\alpha i}.
% y is the vector of gas macroparameters on the current solving step,
% T0 is the characteristic temperature; n0 is the characteristic number 
% density; Prcl is the container for all particles; Coll is the container
% for all collision parameters; setup is the problem setup;
% ind_exc indicates if CO is excited; i_dis is the dissociation model 
% indicator; ind_Arr is the Arrenius constants indicator;
% ind_U is the U parameter indicator in dissociation models;
% ind_Aliat is the Aliat or MT switcher.
% 09.12.2022 by Maksim Melnik.

num_vibr_levels=Prcl.CO.num_vibr_levels(1);
num_COa=num_vibr_levels+1;
num_COA=num_COa+Prcl.CO.num_vibr_levels(2);
num_C=num_COA+Prcl.CO.num_vibr_levels(3);
num_O=num_C+1;
num_C2=num_O+1;
num_Ar=num_C2+1;
num_v=num_Ar+1;
num_T=num_v+1;
R=zeros(num_Ar, 1);
ni_b = y(1:num_vibr_levels);
nCOa_b = y(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
nCOA_b = y(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
nac_b = y(num_C);
nao_b = y(num_O);
nC2_b = y(num_C2);
naAr_b = y(num_Ar);
nm_b = sum(ni_b);
T_b = y(num_T);

    % all zero
R_VT_CO_C2=0;   R_VT_CO_Ar=0;   R_VT_data_COa=0;    R_VT_data_COA=0;
R_VE_CO_toCOa=0;    R_VE_CO_toCOA=0;    R_d_CO_CO=0;    R_d_CO_C2=0;
R_d_CO_Ar=0;    R_d_C2=0;   R_exch_CO_C_C2_O=0;

if i_dis<4
    R_d_CO_CO=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], ...
            nac_b, nao_b, nm_b+sum(nCOa_b)+sum(nCOA_b), Coll.CO_CO, ...
                                 T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
else
 ind_exc_Savelev=3; % 1 - только диссоциирующий атом возбуждён
                    % 2 - ещё + партнёр, 3 - ещё + продукт
    R_d_CO_CO_Savelev=R_diss_Savelev(Prcl.CO, [ni_b; nCOa_b; nCOA_b], ...
        nac_b, nao_b, nm_b+sum(nCOa_b)+sum(nCOA_b), Coll.CO_CO, ...
                            T_b*T0, n0, ind_Arr, ind_U, ind_exc_Savelev);
    R_d_CO_CO=sum(R_d_CO_CO_Savelev, [2, 3]);
end
R_d_CO_C= R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b, ...
        nao_b, nac_b, Coll.CO_C, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
R_d_CO_O= R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b, ...
        nao_b, nao_b, Coll.CO_O, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
if naAr_b>0
 R_d_CO_Ar=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b,...
      nao_b, naAr_b, Coll.CO_Ar, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
end
if setup.C2
 R_d_CO_C2=R_diss_Aliat_onoff(Prcl.CO, [ni_b; nCOa_b; nCOA_b], nac_b,...
      nao_b, nC2_b, Coll.CO_C2, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
end
R_d_data=R_d_CO_CO+R_d_CO_C+R_d_CO_O+R_d_CO_C2+R_d_CO_Ar;
    % COa dissociation
R_d_COa=R_d_data(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
    % COA dissociation
R_d_COA=R_d_data(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
    % C2 dissociation
if setup.C2
 R_d_C2=R_diss_Aliat_onoff(Prcl.C2, nC2_b, nac_b, nac_b, ...
      nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, ...
                        Coll.C2, T_b*T0, n0, ind_Arr, ind_U, ind_Aliat);
 R_exch_CO_C_C2_O=R_exch_CO_C__C2_O(Prcl.C2, Prcl.CO, ...
     [ni_b; nCOa_b; nCOA_b], nac_b, nao_b, nC2_b, Coll.CO_C__C2_O, ...
                                                       T_b*T0, ind_Arr);
end

R_VT_CO_CO=R_VT_old(Prcl.CO, ni_b, Prcl.CO, ...
              nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 1, setup.model_VT);
R_VT_CO_C = ...
    R_VT_old(Prcl.CO, ni_b, Prcl.C,  nac_b,  T_b*T0, 1, setup.model_VT);
R_VT_CO_O = ...
    R_VT_old(Prcl.CO, ni_b, Prcl.O,  nao_b,  T_b*T0, 1, setup.model_VT);
if naAr_b>0
 R_VT_CO_Ar= ...
     R_VT_old(Prcl.CO, ni_b, Prcl.Ar, naAr_b, T_b*T0, 1, setup.model_VT);
end
if setup.C2
 R_VT_CO_C2= ...
     R_VT_old(Prcl.CO, ni_b, Prcl.C2, nC2_b, T_b*T0, 1, setup.model_VT);
end
R_VT_data=R_VT_CO_CO+R_VT_CO_C+R_VT_CO_O+R_VT_CO_C2+R_VT_CO_Ar;

if ind_exc
 R_VT_data_COa=R_VT_old(Prcl.CO, nCOa_b, Prcl.CO, ...         CO, COa, COA
            nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 2, setup.model_VT) ...
  +R_VT_old(Prcl.CO, nCOa_b, Prcl.C, nac_b, T_b*T0, 2, setup.model_VT)...C
  +R_VT_old(Prcl.CO, nCOa_b, Prcl.O, nao_b, T_b*T0, 2, setup.model_VT)...O
  ;
 if naAr_b>0
  R_VT_data_COa=R_VT_data_COa...
   +R_VT_old(Prcl.CO, nCOa_b, Prcl.Ar, naAr_b, T_b*T0, 2, setup.model_VT);
 end
 if setup.C2
  R_VT_data_COa=R_VT_data_COa+...
    R_VT_old(Prcl.CO, nCOa_b, Prcl.C2, nC2_b, T_b*T0, 2, setup.model_VT);
 end
 R_VT_data_COA=R_VT_old(Prcl.CO, nCOA_b, Prcl.CO, ...         CO, COa, COA
            nm_b+sum(nCOa_b)+sum(nCOA_b), T_b*T0, 3, setup.model_VT)...
  +R_VT_old(Prcl.CO, nCOA_b, Prcl.C, nac_b, T_b*T0, 3, setup.model_VT)...C
  +R_VT_old(Prcl.CO, nCOA_b, Prcl.O, nao_b, T_b*T0, 3, setup.model_VT)...O
  ;
 if naAr_b>0
  R_VT_data_COA=R_VT_data_COA...
   +R_VT_old(Prcl.CO, nCOA_b, Prcl.Ar, naAr_b, T_b*T0, 3, setup.model_VT);
 end
 if setup.C2
  R_VT_data_COA=R_VT_data_COA+...
     R_VT_old(Prcl.CO, nCOA_b, Prcl.C2, nC2_b, T_b*T0, 3, setup.model_VT);
 end

    % CO VE
 R_VE_CO_toCOa=R_VE_m(Prcl.CO, ni_b, nCOa_b, ...
    nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, T_b*T0, 2);
 R_VE_CO_toCOA=R_VE_m(Prcl.CO, ni_b, nCOA_b, ...
    nm_b+sum(nCOa_b)+sum(nCOA_b)+nac_b+nao_b+nC2_b+naAr_b, T_b*T0, 3);
end

R(1:num_vibr_levels)=R_VT_data + R_d_data(1:num_vibr_levels)+...
     sum(R_VE_CO_toCOa,2)+sum(R_VE_CO_toCOA,2);
R(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1)= R_d_COa +R_VT_data_COa...
                                                -sum(R_VE_CO_toCOa,1)';
R(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1)= R_d_COA +R_VT_data_COA...
                                                -sum(R_VE_CO_toCOA,1)';
R(num_C)=-sum(R_d_data)-2*R_d_C2+sum(R_exch_CO_C_C2_O);
R(num_O)=-sum(R_d_data)-sum(R_exch_CO_C_C2_O);
end