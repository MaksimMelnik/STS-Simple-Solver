function check_CL_SW(L0, Y, kinetics, ind_ref)
% The function aimed to check if conservation laws are respected in the SW
% problem.
% L0 is the initial values, Y is the evolution of macroparameters, 
% kinetics is the kinetics container, ind_ref is the indicator for
% the reference value (0 -- before SW, 1 -- just behind SW).
% 25.12.2022 Maksim Melnik

    % constants
k=1.380649e-23;             % Boltzmann constant, J/K
rho=0;
for ind=1:kinetics.num_Ps
 rho=rho + sum(Y(:, kinetics.index{ind}), 2) * kinetics.Ps{ind}.mass;
end
rhov=rho.*Y(:, end-1);                                      % rho*v
rhov0=L0(1);
if ind_ref==1
    rhov0=rhov(1);
end
disp(['rho*v max error ' num2str(max(abs((rhov-rhov0)/rhov0)))])
pres=sum(Y(:, 1:end-2), 2)*k.*Y(:, end);                    % p
rhov2p=rho.*Y(:, end-1).^2+pres;                            % rho*v^2+p
rhov2p0=L0(2);
if ind_ref==1
    rhov2p0=rhov2p(1);
end
disp(['rho*v^2+p max error ' num2str(max(abs((rhov2p-rhov2p0)/rhov2p0)))])
En2=0;
for ind=1:kinetics.num_Ps
 if kinetics.Ps{ind}.fr_deg_c>3
  e_i=[];
  for ind_e=1:kinetics.Ps{ind}.num_elex_levels
   e_i=[e_i, kinetics.Ps{ind}.ev_i{ind_e}+kinetics.Ps{ind}.ev_0(ind_e)+...
                    kinetics.Ps{ind}.e_E(ind_e)];
  end
  En2=En2+sum(Y(:, kinetics.index{ind}).*...
                        (e_i+kinetics.Ps{ind}.form_e+2.5*k*Y(:, end)), 2);
 else
  En2=En2+sum(Y(:, kinetics.index{ind}).*...
       (kinetics.Ps{ind}.e_E(1:kinetics.Ps{ind}.num_elex_levels)+...
                            kinetics.Ps{ind}.form_e+1.5*k*Y(:,end)), 2);
 end
end
Ep2=(En2+pres)./rho+0.5*Y(:, end-1).^2;
Ep0=L0(3);
if ind_ref==1
    Ep0=Ep2(1);
end
disp(['E cons max error ' num2str(max(abs((Ep2-Ep0)/Ep0)))])
end