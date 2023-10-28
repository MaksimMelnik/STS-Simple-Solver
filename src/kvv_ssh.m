function kdown = kvv_ssh(t, M1, M2, ind_e1, ind_e2, anharm)
% rate coef-s of VV exchanges, AB - CD
% kdown - array of k_(i->i-1)^(j->j+1)
% t - temperature, K
% Modification for an universal mixture.
% t is the temperature; M1 is the first molecule; M2 is the partner 
% particle; ind_e1 is the number of the electronic state of M1; 
% ind_e2 is the number of the electronic state of M2; 
% anharm is the indicator of anharmonicity.
% 29.12.2022 Maksim Melnik

c = 299792458;                          % speed of light
h = 6.626070041e-34;                    % Plank constant, J*sec
k = 1.380649e-23;                       % Boltzmann constant, J/K
mu = M1.mass*M2.mass/(M1.mass+M2.mass); % reduced mass of A2-A2
m_osc=M1.red_osc_mass;                  % reduced mass of osc-or
om_e=M1.we(ind_e1);                     % \omega_e of M1
om_x_e=M1.wexe(ind_e1);                 % \omega_e x_e of M1
l=M1.num_vibr_levels(ind_e1);           % number of vibr. lvls

r0=(M1.diameter+M2.diameter)/2;

if prod(M1.name==M2.name) && ind_e1 == ind_e2
 if anharm == 1 % number of vibr. levels
    om10 = om_e-2*om_x_e; % lenear oscillation frequency, anh.os.
 elseif anharm == 2
    om10 = om_e; % har.os.
 end

Lmax = l-1; % max level
om0 = 2*pi*c*om10;% circular oscillation frequency, sec^-1
%
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
else    % non-resonance VV-exchange (VV')
 epsm=sqrt(M1.EM*M2.EM);    % LJ depth in K
 f=exp(epsm/t);             % DN
 R1=(1:(M1.num_vibr_levels(ind_e1)-1))';
 R2=1:(M2.num_vibr_levels(ind_e2)-1);
 alpha=(M1.BMbeta+M2.BMbeta)/2;
 Zv1=R1*16.84*alpha^2/(M1.m_mass*(M1.we(ind_e1)-2*M1.wexe(ind_e1))/100);
 Zv2=R2*16.84*alpha^2/(M2.m_mass*(M2.we(ind_e2)-2*M2.wexe(ind_e2))/100);
 mu=M1.m_mass*M2.m_mass/(M1.m_mass+M2.m_mass);
 dOmega=abs(M1.ev_i{ind_e1}(2:end)'+M2.ev_i{ind_e2}(1:end-1) ...
            -M1.ev_i{ind_e1}(1:end-1)'-M2.ev_i{ind_e2}(2:end))/h/c/100;
 thetaLT=0.212*mu*dOmega.^2/alpha^2;
 theta_ast=1.4388*dOmega;
 fun=@(z) exp(-z)./...
  (sinh(2.*sqrt(thetaLT)./theta_ast.*(sqrt(z*t+theta_ast)-sqrt(z*t)))).^2;
 Ztr=1.62*(thetaLT./theta_ast).^2.*integral(fun,0,Inf,'ArrayValued',true);
 Q=1/9*f.*Zv1.*Zv2.*Ztr;
%  kdown=Q;

 sigma=r0*1e10;
 Z=4.571e-12*sqrt(t/mu)*sigma^2;
 kdown=Z.*Q/1e6;
end