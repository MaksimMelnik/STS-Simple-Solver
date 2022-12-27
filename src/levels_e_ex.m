function out=levels_e_ex(M, elvl, ind)
V_H = 6.626070041e-34; V_C = 299792458;
    i=0:M.num_vibr_levels(elvl)-1;
    if nargin==2
       ind=1;
    end
    switch ind
        case 1  % à.î.
            e_i=(M.we(elvl)*(i + 0.5) - M.wexe(elvl)*(i + 0.5).^2 + ...
                    M.weye(elvl)*(i + 0.5).^3)*V_H*V_C;
        case 2  % ã.î.
            e_i=(M.we(elvl)*(i + 0.5))*V_H*V_C;
    end
	out=e_i;
end