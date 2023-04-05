function out=levels_e_ex(M, elvl, ind)
% Vibrational energy calculation for diatomic molecule, accounting
% possibility of electronic excitation.
% M is the molecule, elvl is the electronic level, ind is switcher of
% an anharmonicity.
% Return an array of vibrational energy values, J.
% 10.08.2022
h = 6.626070041e-34;    % kg*m^2/s (J*s)
c = 299792458;          % m/s
    i=0:M.num_vibr_levels(elvl)-1;
    if nargin==2
       ind=1;           % anharmonicity parameter, a.o. by default
    end
    switch ind
        case 1          % anharmonic oscillator
            e_i=(M.we(elvl)*(i + 0.5) - M.wexe(elvl)*(i + 0.5).^2)*h*c;% J
        case 2          % harmonic oscillator
            e_i=(M.we(elvl)*(i + 0.5))*h*c;             % J
        case 3          % advanced anharmonic oscillator
            e_i=(M.we(elvl)*(i + 0.5) - M.wexe(elvl)*(i + 0.5).^2 + ...
                    M.weye(elvl)*(i + 0.5).^3)*h*c;     % J
    end
	out=e_i;
end