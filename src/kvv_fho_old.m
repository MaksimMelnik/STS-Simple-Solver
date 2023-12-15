function kvv = kvv_fho_old(T, M1, M2, ind_e1, ind_e2)
% unfinished
% Rate coefficients of VV exchanges, AB - CD
% kdown - array of k_(i->i-1)^(j->j+1)
% T - temperature, K
% Modification for an universal mixture.
% T is the temperature; M1 is the first molecule; M2 is the partner 
% particle; ind_e1 is the number of the electronic state of M1; 
% ind_e2 is the number of the electronic state of M2.
% 30.01.2023 Maksim Melnik, Olga Kunova

Lmax=M1.num_vibr_levels(ind_e1)-1;
Lmax2=M2.num_vibr_levels(ind_e2)-1;
    % constants
k = 1.380649e-23;         % Boltzmann constant, J/K
c = 29979245800;            % speed of light, cm/s (or m/s should be?)
alpha = 4e10;

gaskin = ((M1.diameter + M2.diameter) / 2 * 1e2) ^ 2;
mu = M1.mass*M2.mass/(M1.mass+M2.mass); % kg
d8pi = sqrt(8*pi);
velt = sqrt(k*T/mu);
freq = gaskin*100*d8pi*velt;

Evibi1 = M1.ev_i{ind_e1}(2:end)'/k;
Evibf1 = M1.ev_i{ind_e1}(1:end-1)'/k;
Evibi2 = M2.ev_i{ind_e2}(1:end-1)/k;
Evibf2 = M2.ev_i{ind_e2}(2:end)/k;

ns1 = 1:Lmax;
ns2 = 1:Lmax2;

    eee1 = abs((Evibi1-Evibf1));
    eee2 = abs((Evibi2-Evibf2));
    Evib = (eee1+eee2)/2;
    
delE = Evibi1-Evibf1+Evibi2-Evibf2;
omeg = 2*pi*Evib*c/1.439;

% ro2 = alpha^2*k*T/(2*omeg^2*mu)/27; % вот это не знаю
ro2 = alpha^2*k*T./(2*omeg.^2*mu)/16;       
 
    % Resonance defect correction (Keck, Carrier, Adamovich)
popc = 4*pi^2 * mu * omeg.^2 / alpha^2 / k;
popc = sqrt(popc/T) .* abs(delE)./Evib / sqrt(8);
popc = (2/3)*popc;
popc = exp(-2/3 * popc);
popc = 0.5*(3-popc).*popc;

    % FHO rates for the exact resonance with the correction
zompl0 = (ns1)' .* ns2;
zompl1 = zompl0.*ro2;
zompl2 = zompl1;
zompl3 = (1+zompl1).^2;
zompl4 = zompl2./zompl3;
kvv = zompl4 .* popc * freq;
    %вот ето не факт, но у јдамовича так
kvv=kvv.*exp(0.5*delE/T)/1e6;
end