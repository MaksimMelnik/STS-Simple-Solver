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
% for .type = "A(T/d_T)^n"
%  - .A
%  - .n
%  - .d_T
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

    % template
template.name = NaN;
template.source = NaN;
template.particles = NaN;
template.type = NaN;
template.direction_forward = true;
template.reverse = false;
template.index = NaN;
template.A   = 0;
template.n   = 0;
template.E   = 0;
template.d_T = 1;

     % zero reaction
zero_r.name = 'zero';
react1 = template;
react1.type = "const";
react1.A   = 0;
react1.n   = 0;
zero_r.data = react1;



     % Zeldovich reaction N2 + O -> NO + N
Zeldovich1.name = 'N2 + O -> NO + N';
Zeldovich1.particles = ["N2", "O", "NO", "N"];
   % Savelev's model [6]
react1 = template;
react1.name = Zeldovich1.name;
react1.source = 'Savelev2018';
react1.particles = Zeldovich1.particles;
react1.type = "Heaviside";
react1.direction_forward = true;
react1.reverse = true;
 % all ground states
react1.index = {{1, "all"}, {1, 1}, {1, "all"}, {1, 1}};
react1.A = @(T) 3e-17^(T < 4000)*1.554e-23^(T >= 4000);
react1.n = @(T) 0^(T < 4000)*1.745^(T >= 4000);
react1.E = 37484 * k;   % J
   % Savelev's model [6], but only the ground vib. state of NO
react2 = react1;
react2.source = 'Savelev2018, NO(1)';
 % only ground NO state included
react2.index = {{1, "all"}, {1, 1}, {1, 1}, {1, 1}};
   % from works by V. Guerra [1]
react3 = react1;
react3.source = 'Guerra95';
react3.type = "const";
react3.reverse = false;
react3.index = {{1, 1+13:N2.num_vibr_levels(1)}, {1, 1}, ...
                                                    {1, "all"}, {1, 1}};
react3.A = 1e-13 / 1e6;
react3.n = 0;
react3.E = 0;
   % from works by V. Guerra [1]
react4 = react1;
react4.source = 'Guerra95_reverse';
react4.type = "ATn";
react4.direction_forward = false;
react4.reverse = false;
% react4.index = {1+1:1+5, 1, 1, 1};
react4.index = {{1, 1+3}, {1, 1}, {1, "all"}, {1, 1}};
react4.A = 1.05e-12 / 1e6;
react4.n = 0.5;
react4.E = 0;
   % from works by V. Guerra [1], data from [4]
react5 = react1;
react5.source = 'Kossyi1992';
% react5.type = "ATn";
react5.type = "Starik_test";
react5.direction_forward = false;
react5.reverse = true;
% react5.reverse = false;
react5.index = {{1, "all"}, {1, 1}, {1, "all"}, {1, 1}};
react5.A = 1.05e-12 / 1e6;
react5.n = 0.5;
react5.E = 0;
keySet = {react1.source, react2.source, react3.source, react4.source, ...
    react5.source};
valueSet = {react1, react2, react3, react4, react5};
Zeldovich1.data = containers.Map(keySet, valueSet);


     % Zeldovich reaction O2 + N -> NO + O
Zeldovich2.name = 'O2 + N -> NO + O';
Zeldovich2.particles = ["O2", "N", "NO", "O"];
   % Savelev's model [6]
react1 = template;
react1.name = Zeldovich2.name;
react1.source = 'Savelev2018';
react1.particles = Zeldovich2.particles;
react1.type = "Heaviside";
react1.direction_forward = true;
react1.reverse = true;
 % all ground states
react1.index = {{1, "all"}, {1, 1}, {1, "all"}, {1, 1}};
react1.A = @(T) 4e-16^(T < 4000)*3.206e-23^(T >= 4000);
react1.n = @(T) (-0.39)^(T < 4000)*1.58^(T >= 4000);
react1.E = 1449 * k;   % J
   % Savelev's model [6], but only the ground vib. state of NO
react2 = react1;
react2.source = 'Savelev2018, NO(1)';
 % only ground NO state included
react2.index = {{1, "all"}, {1, 1}, {1, 1}, {1, 1}};
   % from works by C D Pintassilgo [2] and V Guerra, data from [4]
react3 = react1;
react3.source = 'Kossyi1992_from_Guerra';
react3.type   = "Arrhenius";
react3.reverse = false;
react3.A      = 1.1e-14 / 1e6;
react3.n      = 1;
react3.E      = 3150 * k; % = E/k => E = 3150*k
react3.index = {{1, "all"}, {1, 1}, {1, "all"}, {1, 1}};
   % from [4], depending on T
react4          = react1;
react4.source   = 'Kossyi1992';
react4.type     = "Starik_test_on_T";
react4.reverse  = true;
react4.A        = @(T) (4.5e-12 / 1e6) ^ (T >= 200 && T <= 300) * ...
                                            (1.1e-14 / 1e6) ^ (T > 300);
react4.n        = @(T) T > 300;
react4.E        = @(T) (3220 ^ (T >= 200 && T <= 300) * ...
                                                    3150 ^ (T > 300)) * k;
react4.index = {{1, "all"}, {1, 1}, {1, "all"}, {1, 1}};
keySet = {react1.source, react2.source, react3.source, react4.source};
valueSet = {react1, react2, react3, react4};
Zeldovich2.data = containers.Map(keySet, valueSet);


     % Diffusion coeff of metastable N2 on the wall N2(A) + wall -> N2(X)
N2A_wall_diffusion.name = 'N2(A) + wall -> N2(X) + wall';
N2A_wall_diffusion.particles = ["N2A", "wall", "N2X", "wall"];
   % from works by C D Pintassilgo [2] and V Guerra, data from [3]
react1 = template;
react1.source = 'Levron1978';
react1.type = "A(T/d_T)^n";
react1.A   = 5e18 * 1e2;
react1.d_T = 300;
react1.n   = 0.5;
keySet = {react1.source};
valueSet = {react1};
N2A_wall_diffusion.data = containers.Map(keySet, valueSet);


     % N2(A) + O2 -> N2(X) + O + O
N2A_O2__N2X_O_O.name = 'N2(A) + O2 -> N2(X) + O + O';
N2A_O2__N2X_O_O.particles = ["N2", "O2", "N2", "O", "O"];
   % from works by C D Pintassilgo [2] and V Guerra
react1 = template;
react1.name = N2A_O2__N2X_O_O.name;
react1.particles = N2A_O2__N2X_O_O.particles;
react1.source = 'Pintassilgo2009';
react1.type = "A(T/d_T)^n";
% react1.reverse      = true;
react1.A   = 1.63e-12 / 1e6;
react1.d_T = 300;
react1.n   = 0.55;
react1.index = {{2, "all"}, {1, "all"}, {1, "all"}, {1, 1}, {1, 1}};
keySet = {react1.source};
valueSet = {react1};
N2A_O2__N2X_O_O.data = containers.Map(keySet, valueSet);


     % N2(B) + O2 -> N2(X) + O + O
N2B_O2__N2X_O_O.name = 'N2(B) + O2 -> N2(X) + O + O';
N2B_O2__N2X_O_O.particles = ["N2", "O2", "N2", "O", "O"];
   % from works by C D Pintassilgo [2] and V Guerra, data from [4]
react1 = template;
react1.name = N2B_O2__N2X_O_O.name;
react1.particles = N2B_O2__N2X_O_O.particles;
react1.source = 'Kossyi1992';
react1.type = "const";
react1.A   = 3e-10 / 1e6;% / 1e2;
react1.index = {{3, "all"}, {1, "all"}, {1, "all"}, {1, 1}, {1, 1}};
keySet = {react1.source};
valueSet = {react1};
N2B_O2__N2X_O_O.data = containers.Map(keySet, valueSet);


%      % N2(B) + N2 -> N2(A) + N2
% N2B_N2__N2A_N2.name = 'N2(B) + N2 -> N2(A) + N2';
% N2B_N2__N2A_N2.particles = ["N2", "N2", "N2", "N2"];
%    % from works by C D Pintassilgo [2] and V Guerra [5]
% react1 = template;
% react1.name = N2B_N2__N2A_N2.name;
% react1.particles = N2B_N2__N2A_N2.particles;
% react1.source = 'Guerra1997';
% react1.type = "const";
% react1.A   = 0.95*3e-11 / 1e6;
% react1.index = {{3, "all"}, {1, "all"}, {2, "all"}, {1, "all"}};
% keySet = {react1.source};
% valueSet = {react1};
% N2B_N2__N2A_N2.data = containers.Map(keySet, valueSet);


     % N2+(X) + O2(X) -> O2+(X) + N2
N2pX_O2X__O2pX_N2.name = 'N2+(X) + O2(X) -> O2+(X) + N2';
N2pX_O2X__O2pX_N2.particles = ["N2+", "O2", "O2+", "N2"];
react1 = template;
react1.name = N2pX_O2X__O2pX_N2.name;
react1.particles = N2pX_O2X__O2pX_N2.particles;
react1.source    = 'Kossyi1992';
react1.type      = "A(T/d_T)^n";
react1.index     = {{1, "all"}, {1, "all"}, {1, "all"}, {1, "all"}};
react1.A         = 6e-11 / 1e6;
react1.n         = - 0.5;
react1.d_T       = 300;
keySet           = {react1.source};
valueSet         = {react1};
N2pX_O2X__O2pX_N2.data = containers.Map(keySet, valueSet);


%% e reactions
     % e + N2+(X) -> N(4S) + N(4S)
e_N2pX__N4S_N4S.name = 'e + N2+(X) -> N(4S) + N(4S)';
e_N2pX__N4S_N4S.particles = ["e", "N2+", "N", "N"];
   % from works by C D Pintassilgo [2] and V Guerra
react1                    = template;
react1.name               = e_N2pX__N4S_N4S.name;
react1.particles          = e_N2pX__N4S_N4S.particles;
react1.source             = 'Pintassilgo2009';
react1.type               = "A(T/d_T)^n";
react1.index              = {{1, 1}, {1, "all"}, {1, 1}, {1, 1}};
react1.A                  = 4.8e-7 / 1e6;    % m3 / s
react1.n                  = - 1;
react1.d_T                = 300;
keySet                    = {react1.source};
valueSet                  = {react1};
e_N2pX__N4S_N4S.data      = containers.Map(keySet, valueSet);


     % e + O2+(X) -> O + O
e_O2pX__O_O.name = 'e + O2+(X) -> O + O';
e_O2pX__O_O.particles = ["e", "O2+", "O", "O"];
   % from works by C D Pintassilgo [2] and V Guerra, data from [4]
react1                    = template;
react1.name               = e_O2pX__O_O.name;
react1.particles          = e_O2pX__O_O.particles;
react1.source             = 'Kossyi1992';
react1.type               = "A(T/d_T)^n";
react1.index              = {{1, 1}, {1, "all"}, {1, 1}, {1, 1}};
react1.A                  = 2e-7 / 1e6;    % m3 / s
react1.n                  = - 1;
react1.d_T                = 300;
keySet                    = {react1.source};
valueSet                  = {react1};
e_O2pX__O_O.data          = containers.Map(keySet, valueSet);


%% from LoKI
        % e+N2(X)->e+2N(4S),Excitation
e_N2X__e_N4S_N4S.name = 'e+N2(X)->e+2N(4S),Excitation';
e_N2X__e_N4S_N4S.particles = ["e", "N2", "e", "N", "N"];
   % from works by C D Pintassilgo [2] and V Guerra
react1                    = template;
react1.name               = e_N2X__e_N4S_N4S.name;
react1.particles          = e_N2X__e_N4S_N4S.particles;
react1.source             = 'LoKI-B steady';
react1.type               = "const";
react1.index              = {{1, 1}, {1, "all"}, {1, 1}, {1, 1}, {1, 1}};
react1.A                  = 7.52016817823127e-18;    % m3 / s
keySet                    = {react1.source};
valueSet                  = {react1};
e_N2X__e_N4S_N4S.data     = containers.Map(keySet, valueSet);


        % e+N2(X)->e+e+N2(+,X),Ionization
e_N2X__e_e_N2pX.name = 'e+N2(X)->e+e+N2(+,X),Ionization';
e_N2X__e_e_N2pX.particles = ["e", "N2", "e", "e", "N2+"];
   % from works by C D Pintassilgo [2] and V Guerra
react1                    = template;
react1.name               = e_N2X__e_e_N2pX.name;
react1.particles          = e_N2X__e_e_N2pX.particles;
react1.source             = 'LoKI-B steady';
react1.type               = "const";
react1.index           = {{1, 1}, {1, "all"}, {1, 1}, {1, 1}, {1, "all"}};
react1.A                  = 1.14563748590826e-18;    % m3 / s
keySet                    = {react1.source};
valueSet                  = {react1};
e_N2X__e_e_N2pX.data      = containers.Map(keySet, valueSet);


    % summarizing all reactions in the one container and file
keySet = {zero_r.name, Zeldovich1.name, Zeldovich2.name, ...
 N2A_wall_diffusion.name, N2A_O2__N2X_O_O.name, N2B_O2__N2X_O_O.name, ...
 N2pX_O2X__O2pX_N2.name, e_N2pX__N4S_N4S.name, e_O2pX__O_O.name, ...
 e_N2X__e_N4S_N4S.name, e_N2X__e_e_N2pX.name};
valueSet = {zero_r.data, Zeldovich1.data, Zeldovich2.data, ...
 N2A_wall_diffusion.data, N2A_O2__N2X_O_O.data, N2B_O2__N2X_O_O.data, ...
 N2pX_O2X__O2pX_N2.data, e_N2pX__N4S_N4S.data, e_O2pX__O_O.data, ...
 e_N2X__e_N4S_N4S.data, e_N2X__e_e_N2pX.data};
Reactions = containers.Map(keySet, valueSet);
save reactions.mat Reactions

    % references
% [1] V Guerra et al 1995 J. Phys. D: Appl. Phys. 28 1903
% [2] C D Pintassilgo et al Plasma Sources Sci. Technol. 18 (2009) 025005
% [3] D Levron et al J. Chem. Phys. 69, 2260 (1978); doi: 10.1063/1.436788
% [4] Kossyi I A et al 1992 Plasma Sources Sci. Technol. 1 207
% [5] V Guerra et al 1997 Plasma Sources Sci. Technol. 6 373
% [6] Kustova E et al AIP Conference Proceedings 1959, 060010 (2018)
end