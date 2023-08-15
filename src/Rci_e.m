function [R, Qin] = Rci_e(y, kinetics)
% Universal function for relaxation terms R_{c\alpha i} for processes 
% involving free electrons e-.
% R is the relaxation term Rci_e, Qin is the energy flux.
% y is the vector of gas macroparameters on the current solving step;
% kinetics is the big structure with all kinetics.
% 11.08.2023 by Maksim Melnik.

T = y(end) * kinetics.T0;
R  = zeros(kinetics.num_eq + 1, 1);
Qin = 0;
Te = kinetics.Te;  % K (temporal using way)

% if isKey(kinetics.reactions, 'free_e') % processes involving free electrons
IndexOfMolecules = kinetics.IndexOfMolecules;
e.form_e = 0;
e.num_vibr_levels = 0;
e.e_E             = 0;
e.fr_deg_c        = 3;
IndexOfMolecules("e") = kinetics.num_Ps + 1;
free_e_reactions      = kinetics.reactions("free_e");
index                 = [kinetics.index, kinetics.num_eq + 1];
Ps                    = [kinetics.Ps, e];
for ind_free_e = 1:length(free_e_reactions)
 reaction = free_e_reactions(ind_free_e);
 IOM_M1 = IndexOfMolecules(reaction.particles(1));
 IOM_M2 = IndexOfMolecules(reaction.particles(2));
 IOM_M3 = IndexOfMolecules(reaction.particles(3));
 IOM_M4 = IndexOfMolecules(reaction.particles(4));
 indM1  = index{IOM_M1};
 indM2  = index{IOM_M2};
 indM3  = index{IOM_M3};
 indM4  = index{IOM_M4};
 switch length(reaction.particles)
  case 4
   [R_exch_temp, Q_exch] = R_exch(Ps{IOM_M1}, Ps{IOM_M2}, Ps{IOM_M3}, ...
       Ps{IOM_M4}, y(indM1), y(indM2), y(indM3), y(indM4), Te, reaction);
   R(indM1) = R(indM1) + sum(R_exch_temp, [2 3 4]);
   R(indM2) = R(indM2) + sum(R_exch_temp, [1 3 4])';
   R(indM3) = R(indM3) - reshape(sum(R_exch_temp, [1 2 4]), [], 1);
   R(indM4) = R(indM4) - reshape(sum(R_exch_temp, [1 2 3]), [], 1);
  case 5
   IOM_M5 = IndexOfMolecules(reaction.particles(5));
   indM5  = kinetics.index{IOM_M5};
   [R_exch_temp, Q_exch] = R_exch_23(...
        Ps{IOM_M1}, Ps{IOM_M2}, Ps{IOM_M3}, Ps{IOM_M4}, Ps{IOM_M5}, ...
        y(indM1), y(indM2), y(indM3), y(indM4), y(indM5), T, reaction);
   R(indM1) = R(indM1) + sum(R_exch_temp, [2, 3]);
   R(indM2) = R(indM2) + sum(R_exch_temp, [1, 3])';
   R(indM3) = R(indM3) - reshape(sum(R_exch_temp, [1, 2]), [], 1);
   R(indM4) = R(indM4) - reshape(sum(R_exch_temp, [1 2 3 5]), [], 1);
   R(indM5) = R(indM5) - reshape(sum(R_exch_temp, [1 2 3 4]), [], 1);
  otherwise
       error(["Reactions with current number of particles are not " ...
                                                    "implemented yet"])
 end
 Qin = Qin + Q_exch;
end

Ps_i_O2p = IndexOfMolecules("O2+");
O2p = Ps{Ps_i_O2p};
O2 = Ps{IndexOfMolecules("O2")};

[R_ions_wall_data, Q_ions_wall_data] = ...
                    R_ions_wall(O2p, O2, y(index{Ps_i_O2p}), T, kinetics);
R_ions_wall_data = R_ions_wall_data / kinetics.n0;
Q_ions_wall_data = Q_ions_wall_data / kinetics.n0;
R(index{IndexOfMolecules("O2")}) = R(index{IndexOfMolecules("O2")}) + ...
    R_ions_wall_data;
R(index{IndexOfMolecules("O2+")}) = R(index{IndexOfMolecules("O2+")}) ...
    - R_ions_wall_data;
R(index{IndexOfMolecules("e")}) = R(index{IndexOfMolecules("e")}) ...
    - R_ions_wall_data;
Qin = Qin + Q_ions_wall_data;
end