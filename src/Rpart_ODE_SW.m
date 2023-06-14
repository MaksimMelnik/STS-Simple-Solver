function out=Rpart_ODE_SW(x, y2, kinetics)%#ok<INUSL>
% Right part function for ODE systems for shock waves (SW).
% x is the x coordinate; y is the vector of macroparameters;
% kinetics is the big structure with all kinetics.
% 19.12.2022 Maksim Melnik
    
    % constants
k=1.380648528e-23;
v_DN=y2(end-1);     % dimentionless gas velocity
T_DN=y2(end);

    % relaxation terms
[R, ~] = Rci(y2, kinetics);
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
M2(end, end) = 1.5*sum(y2(1:end-2))+nm;
M2sp=sparse(M2);
R=[R; 0; 0]*kinetics.n0*kinetics.Delta/kinetics.v0;
out=M2sp^(-1)*R;
end