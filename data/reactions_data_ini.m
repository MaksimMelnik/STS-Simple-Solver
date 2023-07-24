function reactions_data_ini
% The file with the data for chemical reactions, included in the present
% program code.
% 21.07.23 Maksim Melnik

% Reaction
% .name 
% . container of different sources
% .particles
% subReaction
% .name
% .particles
% .type = const, ATn, Arrhenius, Heaviside
% .direction_forward = 1 (forward), 0 (backward)
% .reverse = on, off
% .source = name of the model or authors 
% .index
% for .type = "const"
%  - .A
% for .type = "ATn"
%  - .A
%  - .n
% for .type = "Arrhenius"
%  - .A
%  - .n
%  - .E
% for .type = "Heaviside"
%  - .A(T)
%  - .n
%  - .E

load('particles.mat') %#ok<LOAD>
k = 1.380649e-23;         % Boltzmann constant, J/K

     % Zeldovich reaction N2 + O -> NO + N
Zeldovich1.name = 'N2 + O -> NO + N';
Zeldovich1.particles = ["N2", "O", "NO", "N"];
   % from works by O. Kunova
react1.name = Zeldovich1.name;
react1.source = 'Kunova';
react1.particles = Zeldovich1.particles;
react1.type = "Heaviside";
react1.direction_forward = true;
react1.reverse = true;
 % all ground states
react1.index = {1:N2.num_vibr_levels(1), 1, 1:NO.num_vibr_levels(1), 1};
react1.A = @(T) 3e-17^(T < 4000)*1.554e-23^(T >= 4000);
react1.n = @(T) 0^(T < 4000)*1.745^(T >= 4000);
react1.E = 37484 * k;   % J
   % from works by O. Kunova, but only the ground vibrational state of NO
react2 = react1;
react2.source = 'Kunova, NO(1)';
 % only ground NO state included
react2.index = {1:N2.num_vibr_levels(1), 1, 1, 1};
   % from works by V. Guerra [1]
react3 = react1;
react3.source = 'Guerra95';
react3.type = "const";
react3.reverse = false;
react3.A = 1e-13 / 1e6;
react3.n = 0;
react3.E = 0;
react3.index = {1+13:N2.num_vibr_levels(1), 1, 1, 1};
keySet = {react1.source, react2.source, react3.source};
valueSet = {react1, react2, react3};
Zeldovich1.data = containers.Map(keySet, valueSet);

keySet = {Zeldovich1.name};
valueSet = {Zeldovich1.data};
Reactions = containers.Map(keySet, valueSet);

save reactions.mat Reactions

    % references
% [1] V Guerra et al 1995 J. Phys. D: Appl. Phys. 28 1903
end