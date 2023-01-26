function kdown = kvv_ssh(t, M1, M2, ind_e, anharm)
% rate coef-s of VV exchanges, A2 - A2
% kdown - array of k_(i->i-1)^(j->j+1)
% t - temperature, K
% Modification for a universal mixture.
% t is the temperature; M1 is the first molecule; M2 is the parthner 
% particle; ind_e is the number of electronic state of M1; 
% anharm is the indicator of anharmonicity.
% 29.12.2022 Maksim Melnik

c = 299792458;            % speed of light
h = 6.626070041e-34;      % Plank constant, J*sec
k = 1.380649e-23;         % Boltzmann constant, J/K
mu = M1.mass*M2.mass/(M1.mass+M2.mass); % reduced mass of A2-A2
m_osc=M1.red_osc_mass;    % reduced mass of osc-or
om_e=M1.we(ind_e);        % \omega_e of M1
om_x_e=M1.wexe(ind_e);    % \omega_e x_e of M1
l=M1.num_vibr_levels(ind_e);    % number of vibr. lvls

if anharm == 1 % number of vibr. levels
    om10 = om_e-2*om_x_e; % lenear oscillation frequency, anh.os.
elseif anharm == 2
    om10 = om_e; % har.os.
end

Lmax = l-1; % max level
om0 = 2*pi*c*om10;% circular oscillation frequency, sec^-1
%
r0=(M1.diameter+M2.diameter)/2;
a = 17.5 / r0(1); % inverse radius in 1st approximation m^-1
Z = r0(1)^2 / sqrt(mu) * sqrt(8*pi*k*t); % collision frequency, m^3/sec
%
Q10 = 0.5^4 / m_osc * a^2 / om0^2 * 4 * k * t;
k10 = Q10 * Z;

kdown = zeros(Lmax,Lmax);
j_up = 0:Lmax-1;
if anharm == 2
    for i_down = 1:Lmax
	kdown(i_down,:) = i_down * (j_up+1) * k10;
    end
elseif anharm == 1
    % anharmonicity factor for VV transitions
    aA = a * 1e-10; % A^-1
    mu_amu = mu / 1.6605e-27; % amu
    dE = om_x_e*1.4388e-2; % K
    delta = 0.427 / aA * sqrt(mu_amu / t) * dE; % Gorgietz book
    for i_down = 1:Lmax
        kdown(i_down,:) = i_down * (j_up+1) * k10 .* exp(-delta .* ...
            abs(i_down-1-j_up)) .*(1.5-0.5 * exp(-delta .* ...
            abs(i_down-1-j_up))) .* exp((j_up-i_down+1) * h * c * ...
            om_x_e / (k * t));
        if i_down==1
            0.427 / 1e-10 * sqrt(1/ 1.6605e-27 ) * dE; %вот тут не понял
        end
    end
end