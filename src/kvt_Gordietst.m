function kvt = kvt_Gordietst(T, M1, ind_e)
% VT rate coefficient function for N2-O collision presented in [1].
% T is the gas temperature; M1 is N2 molecule; 
% ind_e is the number of electronic state of M1.
% 17.07.2023 Maksim Melnik
% [1] B Gordiets et al. Plasma Sources Sci. Technol. 2 (1993) 158-163

P10 = 2.3e-13*exp(-1280/T) + 2.7e-11*exp(-10840/T); %1->0 rate coefficient
kvt = P10 * (1:M1.num_vibr_levels(ind_e)-1);
kvt = kvt / 1e6;
end