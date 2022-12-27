function out=levels_e(num_vibr_levels, ind)
V_H = 6.626070041e-34; V_C = 299792458;
V_ome = 216981; V_omexe = 1329; V_omeye=1.0511;
    i=0:num_vibr_levels-1;
    if nargin==1
       ind=1;
    end
    switch ind
        case 1
            e_i=(V_ome*(i + 0.5) - V_omexe*(i + 0.5).^2 + ...   %à.î.
                    V_omeye*(i + 0.5).^3)*V_H*V_C;
        case 2
            e_i=(V_ome*(i + 0.5))*V_H*V_C;                      %ã.î.
    end
	out=e_i;
end