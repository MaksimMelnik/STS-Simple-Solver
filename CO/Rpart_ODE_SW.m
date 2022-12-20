function out=Rpart_ODE_SW(x, y, Prcl, ind_exc, setup, kinetics)%#ok<INUSL>
% Right part function for ODE systems for shock waves (SW).
% x is the x coordinate; y is the vector of macroparameters;
% Prcl is the particles container;
% ind_exc is the indicator of the electronic excitation 
% presence; setup is the different parameters for kinetic scheme container;
% kinetics is the big structure with all kinetics.
% 19.12.2022 Maksim Melnik
    
    % constants
k=1.380648528e-23;

% num_vibr_levels=Prcl.CO.num_vibr_levels(1);
% num_COa=num_vibr_levels+1;
% num_COA=num_COa+Prcl.CO.num_vibr_levels(2);
% num_C=num_COA+Prcl.CO.num_vibr_levels(3);
% num_O=num_C+1;
% num_C2=num_O+1;
% num_Ar=num_C2+1;
% ni_b = y(1:num_vibr_levels);
% nCOa_b = y(num_COa:num_COa+Prcl.CO.num_vibr_levels(2)-1);
% nCOA_b = y(num_COA:num_COA+Prcl.CO.num_vibr_levels(3)-1);
% nac_b = y(num_C);
% nao_b = y(num_O);
% nC2_b = y(num_C2);
% naAr_b = y(num_Ar);
% 
%     % call of relaxation terms Rci
% y2=ni_b;
% if ind_exc
%     y2=[y2; nCOa_b; nCOA_b];
% end
% y2=[y2; nac_b; nao_b];
% if setup.C2
%     y2=[y2; nC2_b];
% end
% if setup.f<1
%     y2=[y2; naAr_b];
% end
% y2=[y2; y(end-1:end)];
y2=y;

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
M2(end-1, end-1) = rho*v_DN*kinetics.v0^2/k/kinetics.T0;
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
  M2(end, kinetics.index{ind}) = ...
                (e_i+ kinetics.Ps{ind}.form_e)/k/kinetics.T0 +2.5*T_DN;
  M2(end, end-1)=M2(end, end-1)+...
            (e_i +kinetics.Ps{ind}.form_e)/k/kinetics.T0*...
            y2(kinetics.index{ind})+3.5*sum(y2(kinetics.index{ind}))*T_DN;
  nm=nm+sum(y2(kinetics.index{ind}));
 else
  M2(end, kinetics.index{ind}) = (kinetics.Ps{ind}.form_e ...
      + kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels))...
                                                /k/kinetics.T0 +1.5*T_DN;
  M2(end, end-1)=M2(end, end-1)+2.5*sum(y2(kinetics.index{ind}))*T_DN+...
      (kinetics.Ps{ind}.e_E+kinetics.Ps{ind}.form_e)/k/kinetics.T0*...
      y2(kinetics.index{ind});
 end
end
M2(end, end-1)=M2(end, end-1)/v_DN;
M2(end, end)=1.5*sum(y2(1:end-2))+nm;
M2sp=sparse(M2);
R2=[R2; 0; 0]*kinetics.n0*kinetics.Delta/kinetics.v0;
out_temp=M2sp^(-1)*R2;

if ind_exc==0
    out=[out_temp(1:68); zeros(120-68, 1)];
else
    out=out_temp(1:120);
end
out=[out; out_temp(kinetics.index{2}); out_temp(kinetics.index{3})];
index_last=kinetics.index{3}+1;
if setup.C2
    out=[out; out_temp(index_last)];
    index_last=index_last+1;
else
    out=[out; 0];
end
if setup.f<1
    out=[out; out_temp(index_last)];
else
    out=[out; 0];
end
out=[out; out_temp(end-1:end)];

out=out_temp;

end