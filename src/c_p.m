function c_p = c_p(M, T)
% Function for calculation molar heat capacity c_p (J/mol/K).
% M is the molecule under consideration, T is the gas temperature.
% Data from [1].
% 04.07.2023 by Maksim Melnik.
% [1] C D Pintassilgo et al Plasma Sources Sci. Technol. 23 (2014) 025006.

switch M.name
 case "N2"
  c_p = 29.1 + 2494.2/553.4/sqrt(pi/2) * exp(-2 * ((T-1047.4)/553.4).^2);
 case "O2"
  c_p = 28.2 + 6456.2/788.3/sqrt(pi/2) * exp(-2 * ((T-1006.9)/788.3).^2);
 otherwise
  error(['No data for the heat capasity for the molecule ' M.name '.'])
end
end