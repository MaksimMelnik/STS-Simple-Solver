function k_down = kvt_ssh(t, M1, M2, ind_e, anharm)
% function for vib. exch. rates, SHH model, Stupochenko
% t - temperature, K
% sw_o - oscillator switcher
% k_down - array of the rate coef-s of VT exchanges A2(i+1)+M = A2(i)+M
% i+1 - vib. level before collision
% i - vib. level after collision
% 05.12.2022 modification by Maksim Melnik
% t is the temperature (K), M1 is the first molecule, M2 is the second 
% partner-particle ind_e is the electronic state of M1, anharm is 
% the anharmonicity indicator.

% constants
k = 1.380649e-23;               % Boltzmann constant, J/K
c = 299792458;                  % speed of light (m/s)
h = 6.62607015e-34;             % Plank constant, J*sec
h_bar = h/(2*pi);
om_e=M1.we(ind_e);              % \omega_e of M1
om_x_e=M1.wexe(ind_e);          % \omega_e x_e of M1
l=M1.num_vibr_levels(ind_e);    % number of vibr. lvls

om10 = om_e;                    % lenear oscillation frequency
if anharm                       % if anharmonic 
    om10=om10-2*om_x_e;
end
Lmax = l-1; % number of maximum level
i = (0:Lmax-1);
om0 = 2*pi*c*om10;% circular oscillation frequency, sec^-1

% reduced masses, kg
mu=M1.mass*M2.mass/(M1.mass+M2.mass);
% mu = m(1)*m./(m(1)+m);

% inverse radius in 1st approximation m^-1
r0=(M1.diameter+M2.diameter)/2;
a = 17.5./r0;
% potential depth, Lennard-Jones
em=sqrt(M1.EM*M2.EM*M1.diameter^6*M2.diameter^6)/r0^6;

% collision frequency, m^3/sec
Z = r0.^2./sqrt(mu)*sqrt(8*pi*k*t);

% firstly, find rates of transition from 1st lev. to 2nd lev.
chi = (pi^2*om0^2/2/k*(mu./a.^2)*t^(-1)).^(1/3);% dim-less
%
r = r0 .* (0.5 * sqrt(1+(chi*t)./em) + 0.5).^(-1/6);
% steric factor
Z_0 = (a .* M1.r_e).^2 .* exp(-((3*a .* M1.r_e.^2)./(8*r)));

Const = 1.294./Z_0*4*pi^2*om0/h_bar*sqrt(4*pi/3);
%
r_c_r_02 = (0.5*sqrt(1+(chi*t)./em)+0.5).^(-1/3);
%
P10 = Const.*r_c_r_02.*(1+1.1*em/t).^(-1).*...
                mu./a.^2.*chi.^0.5.*exp(-3*chi+h_bar*om0/(2*k*t)+em/t);
k10 = Z.*P10; % m^3/sec

% secondly, find rates coef-s for transition i+1 -> i
if anharm   % in case of anharmonic oscillator
    % anharmonicity factor for VT transitions (Gordietz)
    aA = a * 1e-10; % A^-1
    % reduced masses
    mu_amu = M1.mass*M2.mass./(M1.mass+M2.mass)/1.6605e-27; % amu
    % adiabatic factor for transition i+1 -> i
    % E(i+1) - E(i)
    diffE = (om_e-2*(i+1)*om_x_e)*h*c; % J
    %
    E1 = (om_e-2*om_x_e)*1.4388e-2; % K
    dE = om_x_e*1.4388e-2; % K
    %
    gamma_n = pi ./ a / h_bar .* sqrt(mu/(2*k*t));
    gamma_n = gamma_n * diffE;
    %
    gamma_0 = 0.32 ./ aA .* sqrt(mu_amu/t)*E1;
    %
    delta= (4/3*gamma_0(1) * dE/E1).^(gamma_n < 20).*...
                        (4*(gamma_0(1))^(2/3) * dE/E1).^(gamma_n >= 20);
    %
    k_down = (i+1) * k10 .* exp(i .* delta) .* exp(-i*h*c*om_x_e / (k*t));
else    % harmonic oscillator case
    k_down = k10*(i+1);
end
end
