function out=Rpart_Tiago_imitation(t, y, Prcl, Coll, tau, xspan, ...
    T0, n0, is_e, EdN, fileData, fileKrate)
% правые части для системы ДУ для нульмерной релаксации в кислороде с 
% электронным возбуждением. Хотим повторить результаты Тиаго, поэтому
% делаем закос под его код.
% 26.11.2020
k = 1.3807e-23;  % постоянная Больцмана (Дж/К)
N_a = 6.0221e23; % постоянная Авогадро (1/моль)

O2=Prcl.O2; O=Prcl.O;
num_O2=O2.num_vibr_levels(1);
num_O2a=num_O2+O2.num_vibr_levels(2);
num_O2b=num_O2a+O2.num_vibr_levels(3);
num_O2p=num_O2b+1;
num_O=num_O2p+1;
num_O1D=num_O+1;
num_Op=num_O1D+1;
num_Om=num_Op+1;
num_O3=num_Om+1;
num_O3exc=num_O3+1;
num_e=num_O3exc+1;
num_T=num_e+1;
R=zeros(num_T,1);
nO2_b=y(1:num_O2);
nm=sum(y(1:num_O2b));%sum(nO2_b);
nO_b=y(num_O);
ne_b=y(num_e);
T=y(num_T)*T0;
nAll_b=sum(y(1:num_O3exc));
if nAll_b~=nm+y(num_O2p)+nO_b+y(num_O1D)+y(num_Op)+y(num_Om)+y(num_O3)+...
        y(num_O3exc)
    disp('ALERT')
end

R_VT_data=zeros(num_O2,1);
R_VV_data=zeros(num_O2,1);
R_diss_data=zeros(num_O2b,1);
% R_eD_data=zeros(num_O2+1,1);
R_eV_data=zeros(25,1);
R_rec_O2X_data=0;        R_rec_O2a_data=0;        R_rec_O2b_data=0;
R_rec_3O3P__O2X_O3P_data=0;                       R_ex_O2_O2a_data=0;
R_deex_O2b_O2_byO3p=0;   R_deex_O2b_O2a_byO3p=0;  R_ex_O2a_O2b_data=0;
R_mix_1DX_3Pb_data=0;    R_deex_1D3PbyO2_data=0;  R_ex_O3P_O1D_data=0;
R_d_O2X_3P1D_data=0;     R_d_O2a_3P1D_data=0;     R_d_O2b_O3P_O1D_data=0;
R_d_O3X__O3P_O2X_data=0; R_ion_O2X__e_O2p_data=0; R_ion_O2a__e_O2p_data=0;
R_diss_e_O2p__2O3P_data=0;              R_diss_e_O2p__O3P_O1D_data=0;
R_ion_O3P_Op_data=0; 
R_mix_O2X__e_O3P_Op_data=0;             R_mix_O2a__e_O3P_Op_data=0;
R_mix_e_O2X__Om_O3P_data=0;             R_mix_e_O2a__Om_O3P_data=0;
R_deion_e_Om__2e_O3P_data=0;            R_mix_2O2a_O2XO2b_data=0;
R_deex_O2a_O3P_O2X_O3P_data=0;          R_deex_O2a_O2__O2X_O2_data=0;
R_rec_O3P_O2X_O3P__O3X_O3P_data=0;      R_rec_O3P_O2X_O2__O3X_O2_data=0;
R_d_O2b_O3X__2O2X_O3P_data=0;           R_d_O2a_O3X__2O2X_O3P_data=0;
R_exch_O3P_O3X__2O2X_data=0;            R_exch_O3P_O3X__O2a_O2X_data=0;
R_exch_O3P_O3X__O2b_O2X_data=0;         R_rec_O3P_O2X_O3X__2O3X_data=0;
R_exch_O1D_O3X__2O2X_data=0;            R_exch_O1D_O3X__O2X_2O3P_data=0;
R_exch_Op_O3X__O2p_O2X_data=0;          R_exch_Op_O2X__O2p_O3P_data=0;
R_exch_Op_O2a__O2p_O3P_data=0;          R_mix_Om_O2a__e_O3X_data=0;
R_mix_Om_O3P__e_O2X_data=0;             R_mix_Om_O2X__e_O3X_data=0;
R_mix_Om_O2b__e_O3P_O2X_data=0;         R_deion_Op_Om__2O3P_data=0;
R_deion_O2p_Om__O2X_O3P_data=0;
R_rec_O3P_2O2X__O3exc_O2X_data=0;       R_deex_O3exc_O3P__O3X_O3P_data=0;
R_deex_O3exc_O2X__O3X_O2X_data=0;       R_mix_O2a_O3exc__2O2X_O3P_data=0;
R_exch_O3P_O3exc__2O2X_data=0;

Qin=0;

%     R_VT_data=R_VT(O2, T, nO2_b, nO_b);
%     R_VV_data=R_VV(O2, O2, T, nO2_b, nO2_b);
    % O2(X,a,b)+O2(X,a,b)->2O(3P)+O2(same)
R_diss_data=R_diss_data + R_diss_Aliat(... 
    Prcl.O2, y(1:num_O2b), nO_b, nO_b, nm, Coll.O2_O2, T,n0);
    % O2(X,a,b)+O(3P)->2O(3P)+O(3P)
R_diss_data=R_diss_data + R_diss_Aliat(Prcl.O2, ...
    y(1:num_O2b), nO_b, nO_b, sum(y(num_O:num_Om)), Coll.O2_O, T, n0);

% if is_e~=0
%     if t*tau>1e-3
%         run_LoKI(T, nO2_b, y(num_O2a), y(num_O2b), ...
%             ne_b, nO_b, y(num_O1D), n0, T_e, EdN, t*tau);
%     end
% end
% if mod(round(cputime), 8)==1
%     if is_e~=0
%         run_LoKI(T, nO2_b, y(num_O2a), y(num_O2b), ...
%             ne_b, nO_b, y(num_O1D), n0, EdN, t*tau);
%     end
% end
if is_e~=0
        % пока нет третьего возбужденного состояния
%     run_LoKI(T, nO2_b, y(num_O2a), y(num_O2b), ...
%         ne_b, nO_b, y(num_O1D),...
%         n0, T_e, EdN, fileKrate);
    R_eV_data=R_eV(nO2_b, y(num_O2a), y(num_O2b), ...
        ne_b, nO_b, y(num_O1D), y(num_Om),...
        y(num_O3), EdN, fileKrate, t*tau, n0);
end
    % recombination 2O(3P)+O2(X)->O2(X)+O2(X)     % можно nm вместо nO2_b?
R_rec_O2X_data=0.5*3.81e-30/T*exp(-170/T)/1e12*y(num_O)*y(num_O)*n0*nO2_b;
    % recombination 2O(3P)+O2(X)->O2(a)+O2(X)     % можно nm вместо nO2_b?
R_rec_O2a_data=0.33*3.81e-30/T*exp(-170/T)/1e12*nO_b*nO_b*nO2_b*n0;
    % recombination 2O(3P)+O2(X)->O2(b)+O2(X)     % можно nm вместо nO2_b?
R_rec_O2b_data=0.17*3.81e-30/T*exp(-170/T)/1e12*nO_b*nO_b*nO2_b*n0;
    % recombination 2O(3P)+O(3P)->O2(X)+O(3P)
R_rec_3O3P__O2X_O3P_data=...
                    3.6e-32*(1/T)^0.63/1e12*y(num_O)*y(num_O)*y(num_O)*n0;
    % R_eV_data(1-5) - R_eV_i, i от 0 до 4
k_eD_eq=R_eV_data(6);
% k_el=R_eV_data(7);
    % Te
Te=R_eV_data(end-1)/k*1.6021766208e-19;
    % excitation e+O2(X)<->e+O2(a)
R_ex_O2_O2a_data=R_eV_data(7);
    % excitation e+O2(X)<->e+O2(b)
R_ex_O2_O2b_data=R_eV_data(8);
    % excitation e+O2(a)<->e+O2(b)
R_ex_O2a_O2b_data=R_eV_data(9);
    % excitation e + O(3P) <-> e + O(1D)
R_ex_O3P_O1D_data=R_eV_data(10);
    % diss e + O2(X) -> e + O(3P) + O(1D)
R_d_O2X_3P1D_data=R_eV_data(11);
    % diss e + O2(a1Dg) -> e + O(3P) + O(1D)
R_d_O2a_3P1D_data=R_eV_data(12);
    % diss e + O2(alDg) -> e + 2O(3P)
k_eD_O2a_2O3P=R_eV_data(13);
    % diss e + O2(b) -> e + 2O(3P)
k_eD_O2b_2O3P=R_eV_data(14);
    % diss e + O2(b) -> e + O(3P) + O(1D)
R_d_O2b_O3P_O1D_data=R_eV_data(15);
    % diss e + O3(X) -> e + O(3P) + O2(X)
R_d_O3X__O3P_O2X_data=R_eV_data(16);
    % ion e + O2(X) -> 2e + O2(+,X)
R_ion_O2X__e_O2p_data=R_eV_data(17);        % not affect on ne
    % ion e + O2(a) -> 2e + O2(+,X)
R_ion_O2a__e_O2p_data=R_eV_data(18);        % not affect on ne
    % diss e + O2(+,X) -> 2O(3P)            % not affect on ne
R_diss_e_O2p__2O3P_data=2e-7*(300/Te)/1e6*ne_b*y(num_O2p);
    % diss e + O2(+,X) -> O(3P) + O(1D)     % not affect on ne
R_diss_e_O2p__O3P_O1D_data=1.95e-7*(300/Te)^0.7/1e6*ne_b*y(num_O2p);
    % e + O(3P) -> 2e + O(+,gnd)
R_ion_O3P_Op_data=R_eV_data(19);            % not affect on ne
    % e + O2(X) -> 2e + O(3P) + O(+,gnd)
R_mix_O2X__e_O3P_Op_data=R_eV_data(20);     % not affect on ne
    % e + O2(a1Dg) -> 2e + O(3P) + O(+,gnd)
R_mix_O2a__e_O3P_Op_data=R_eV_data(21);     % not affect on ne
%     % e + O2(X) -> O(-,gnd) + O(3P)       % not affect on ne
R_mix_e_O2X__Om_O3P_data=R_eV_data(22);
%     % e + O2(a1Dg) -> O(-,gnd) + O(3P)
R_mix_e_O2a__Om_O3P_data=R_eV_data(23);
%     % e + O(-,gnd) -> 2e + O(3P)
R_deion_e_Om__2e_O3P_data=R_eV_data(24);
Q_el=R_eV_data(end);

    % deexcitation O2(b)+O(3P)->O2(X)+O(3P)
R_deex_O2b_O2_byO3p =4e-14/1e6*y(num_O2b)*nO_b;
    % deexcitation O2(b)+O(3P)->O2(a)+O(3P)
R_deex_O2b_O2a_byO3p=4e-14/1e6*y(num_O2b)*nO_b;
    % deexcitation O(1D)+O2(X)->O(3P)+O2(X)
R_deex_1D3PbyO2_data=7e-12*exp(67/T)/1e6*y(num_O1D)*y(num_O2);
    % deexcitation O2(a) + O(3P) -> O2(X) + O(3P)
R_deex_O2a_O3P_O2X_O3P_data=7e-17/1e6*y(num_O2a)*y(num_O);
    % deexcitation O2(a) + O2(X) -> O2(X) + O2(X) % можно nm вместо nO2_b?
R_deex_O2a_O2__O2X_O2_data=2.2e-18*(T/300)^0.8/1e6*y(num_O2a)*nO2_b;
    % deexcitation O(1D) + O(3P) -> O(3P) + O(3P)
R_deex_O1D_O3P__2O3P_data=8e-12/1e6*y(num_O1D)*y(num_O);

    % exchange O(1D) + O2(X) -> O(3P) + O2(b1Sg+)
R_mix_1DX_3Pb_data=2.56e-11*exp(67/T)/1e6*y(num_O1D)*sum(nO2_b);
    % exchange O(1D) + O2(X) -> O(3P) + O2(a1Dg)
R_mix_1DX_3Pa_data=1e-12/1e6*y(num_O1D)*y(num_O2);
    % exchange 2O2(a1Dg) -> O2(b1Sg+) + O2(X)
R_mix_2O2a_O2XO2b_data=...
               1.81e-18*exp(700/T)*(T/300)^3.8/1e6*y(num_O2a)*y(num_O2a);

    % O3
    % O(3P) + O2(X) + O(3P) -> O3(X) + O(3P)
R_rec_O3P_O2X_O3P__O3X_O3P_data=...
                            2.1e-34*exp(345/T)/1e12*nO2_b*nO_b*nO_b*n0;
    % O(3P) + O2(X) + O2(X) -> O3(X) + O2(X)      % можно nm вместо nO2_b?
R_rec_O3P_O2X_O2__O3X_O2_data=...
                        0.33*6.4e-35*exp(663/T)/1e12*nO_b*nO2_b*nO2_b*n0;
    % O2(b1Sg+) + O3(X) -> 2O2(X) + O(3P)
R_d_O2b_O3X__2O2X_O3P_data=1.5e-11/1e6*y(num_O2b)*y(num_O3);
    % O2(a1Dg) + O3(X) -> 2O2(X) + O(3P)
R_d_O2a_O3X__2O2X_O3P_data=...
                        5.2e-11*exp(-2840/T)/1e6*y(num_O2a)*y(num_O3);
    % O(3P) + O3(X) -> 2O2(X)
R_exch_O3P_O3X__2O2X_data=...
                        0.5*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
    % O(3P) + O3(X) -> O2(a1Dg) + O2(X)
R_exch_O3P_O3X__O2a_O2X_data=...
                        0.33*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
    % O(3P) + O3(X) -> O2(b1Sg+) + O2(X)
R_exch_O3P_O3X__O2b_O2X_data=...
                        0.17*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
	% O(3P) + O2(X) + O3(X) -> O3(X) + O3(X) (low)
R_rec_O3P_O2X_O3X__2O3X_data=...
                1.66e-34*exp(T/300)/1e12*y(num_O)*y(num_O2)*y(num_O3)*n0;
    % O(1D) + O3(X) -> 2O2(X) (low)
R_exch_O1D_O3X__2O2X_data=1.2e-10/1e6*y(num_O1D)*y(num_O3);
    % O(1D) + O3(X) -> O2(X) + 2O(3P) (low)
R_exch_O1D_O3X__O2X_2O3P_data=1.2e-10/1e6*y(num_O1D)*y(num_O3);

    % O+
    % exch O(+,gnd) + O3(X) -> O2(+,X) + O2(X)      % not affect on ne
R_exch_Op_O3X__O2p_O2X_data=1e-10/1e6*y(num_Op)*y(num_O3);
    % exch O(+,gnd) + O2(X) -> O2(+,X) + O(3P)      % not affect on ne
R_exch_Op_O2X__O2p_O3P_data=2e-11*(300/T)^0.5/1e6*y(num_Op)*y(num_O2);
    % exch O(+,gnd) + O2(a1Dg) -> O2(+,X) + O(3P)   % not affect on ne
R_exch_Op_O2a__O2p_O3P_data=2e-11*(300/T)^0.5/1e6*y(num_Op)*y(num_O2a);

    % O-
    % O(-,gnd) + O2(a1Dg) -> e + O3(X)
R_mix_Om_O2a__e_O3X_data=0.75*1.9e-10/1e6*y(num_Om)*y(num_O2a);
    % O(-,gnd) + O(3P) -> e + O2(X)  тут krate Тиаго, а не из статьи
R_mix_Om_O3P__e_O2X_data=1.9e-10/1e6*y(num_Om)*y(num_O);
    % O(-,gnd) + O2(X) -> e + O3(X)
R_mix_Om_O2X__e_O3X_data=1e-12/1e6*y(num_Om)*y(num_O2);
    % O(-,gnd) + O2(b1Sg+) -> e + O(3P) + O2(X)
R_mix_Om_O2b__e_O3P_O2X_data=6.9e-10/1e6*y(num_Om)*y(num_O2b);
    % O(+,gnd) + O(-,gnd) -> 2O(3P)
R_deion_Op_Om__2O3P_data=2.8e-7/1e6*y(num_Op)*y(num_Om);
    % O2(+,X) + O(-,gnd) -> O2(X) + O(3P)
R_deion_O2p_Om__O2X_O3P_data=9.6e-8*(300/T)^0.5/1e6*y(num_O2p)*y(num_Om);

    % O3(exc)
    % O(3P) + O2(X) + O2(X) -> O3(exc) + O2(X)
R_rec_O3P_2O2X__O3exc_O2X_data=...
            0.67*6.4e-35*exp(663/T)/1e12*y(num_O)*y(num_O2)*y(num_O2)*n0;
    % O3(exc) + O(3P) -> O3(X) + O(3P)
R_deex_O3exc_O3P__O3X_O3P_data=2e-13/1e6*y(num_O3exc)*y(num_O);
    % O3(exc) + O2(X) -> O3(X) + O2(X)
R_deex_O3exc_O2X__O3X_O2X_data=3e-15/1e6*y(num_O3exc)*y(num_O2);
    % O2(a1Dg) + O3(exc) -> 2O2(X) + O(3P)
R_mix_O2a_O3exc__2O2X_O3P_data=...
                26e-11*exp(-(2840-1553)/T)/1e6*y(num_O2a)*y(num_O3exc);
    % O(3P) + O3(exc) -> 2O2(X)
R_exch_O3P_O3exc__2O2X_data=...
                    8e-12*exp(-(2060-1553)/T)/1e6*y(num_O)*y(num_O3exc);
                    
R_eD=R_diss_e_Aliat(Prcl.O2, ...
    y(1:num_O2b), nO_b, nO_b, ne_b, ...
    [k_eD_eq k_eD_O2a_2O3P k_eD_O2b_2O3P], T, n0);

R(1:num_O2)=R_VT_data+R_VV_data;
R(1:num_O2b)=R(1:num_O2b)+R_diss_data+R_eD;                        % O2
    % aditional processes
% если введу уровни для O2, переписать эту строку
R(num_O2)=R(num_O2)-R_ex_O2_O2a_data-R_ex_O2_O2b_data+...          % O2(X)
    R_deex_O2b_O2_byO3p-R_mix_1DX_3Pb_data-R_d_O2X_3P1D_data...
    +R_mix_2O2a_O2XO2b_data-R_rec_O3P_O2X_O3P__O3X_O3P_data...
    -R_rec_O3P_O2X_O2__O3X_O2_data+2*R_d_O2b_O3X__2O2X_O3P_data...
    +2*R_d_O2a_O3X__2O2X_O3P_data+2*R_exch_O3P_O3X__2O2X_data...
    +R_exch_O3P_O3X__O2a_O2X_data+R_exch_O3P_O3X__O2b_O2X_data...
    +R_deex_O2a_O3P_O2X_O3P_data+R_rec_O2X_data...
    +R_rec_3O3P__O2X_O3P_data+R_deex_O2a_O2__O2X_O2_data...
    -R_mix_1DX_3Pa_data+R_d_O3X__O3P_O2X_data...
    -R_rec_O3P_O2X_O3X__2O3X_data+2*R_exch_O1D_O3X__2O2X_data...
    +R_exch_O1D_O3X__O2X_2O3P_data-R_ion_O2X__e_O2p_data...
    +R_exch_Op_O3X__O2p_O2X_data-R_exch_Op_O2X__O2p_O3P_data...
    -R_mix_O2X__e_O3P_Op_data+R_mix_Om_O3P__e_O2X_data...
    -R_mix_Om_O2X__e_O3X_data+R_mix_Om_O2b__e_O3P_O2X_data...
    +R_deion_O2p_Om__O2X_O3P_data-R_mix_e_O2X__Om_O3P_data...
    -R_rec_O3P_2O2X__O3exc_O2X_data+2*R_mix_O2a_O3exc__2O2X_O3P_data...
    +2*R_exch_O3P_O3exc__2O2X_data;
R(num_O2a)=R(num_O2a)+R_rec_O2a_data+R_ex_O2_O2a_data...           % O2a
    +R_deex_O2b_O2a_byO3p-R_ex_O2a_O2b_data-R_d_O2a_3P1D_data...
    -2*R_mix_2O2a_O2XO2b_data-R_d_O2a_O3X__2O2X_O3P_data...
    +R_exch_O3P_O3X__O2a_O2X_data-R_deex_O2a_O3P_O2X_O3P_data...
    -R_deex_O2a_O2__O2X_O2_data+R_mix_1DX_3Pa_data...
    -R_ion_O2a__e_O2p_data-R_exch_Op_O2a__O2p_O3P_data...
    -R_mix_O2a__e_O3P_Op_data-R_mix_Om_O2a__e_O3X_data...
    -R_mix_e_O2a__Om_O3P_data-R_mix_O2a_O3exc__2O2X_O3P_data;
% if t*tau>1e1
%     disp('R')
%     disp(R_rec_O2a_data*n0*n0)
%     disp(R_ex_O2_O2a_data*n0*n0)
%     disp(R_deex_O2b_O2a_byO3p*n0*n0)
%     disp(R_ex_O2a_O2b_data*n0*n0)
%     disp(R_d_O2a_3P1D_data*n0*n0)
%     disp(R_mix_2O2a_O2XO2b_data*n0*n0)
%     disp(R_d_O2a_O3X__2O2X_O3P_data*n0*n0)
%     disp(R_exch_O3P_O3X__O2a_O2X_data*n0*n0)
%     disp(R_deex_O2a_O3P_O2X_O3P_data*n0*n0)
%     disp(R_deex_O2a_O2__O2X_O2_data*n0*n0)
%     disp(R_mix_1DX_3Pa_data*n0*n0)
%     disp(R_ion_O2a__e_O2p_data*n0*n0)
%     disp(R_exch_Op_O2a__O2p_O3P_data*n0*n0)
%     disp(R_mix_O2a__e_O3P_Op_data*n0*n0)
%     disp(R_mix_Om_O2a__e_O3X_data*n0*n0)
%     disp(R_mix_e_O2a__Om_O3P_data*n0*n0)
%     disp(R_mix_O2a_O3exc__2O2X_O3P_data*n0*n0)
% end
R(num_O2b)=R(num_O2b)+R_rec_O2b_data+R_ex_O2_O2b_data...           % O2b
    -R_deex_O2b_O2_byO3p-R_deex_O2b_O2a_byO3p+R_ex_O2a_O2b_data...
    +R_mix_1DX_3Pb_data+R_mix_2O2a_O2XO2b_data...
    -R_d_O2b_O3X__2O2X_O3P_data+R_exch_O3P_O3X__O2b_O2X_data...
    -R_d_O2b_O3P_O1D_data-R_mix_Om_O2b__e_O3P_O2X_data;
R(num_O2p)=R_ion_O2X__e_O2p_data+R_ion_O2a__e_O2p_data...
    -R_diss_e_O2p__2O3P_data-R_diss_e_O2p__O3P_O1D_data...
    +R_exch_Op_O3X__O2p_O2X_data+R_exch_Op_O2X__O2p_O3P_data...
    +R_exch_Op_O2a__O2p_O3P_data-R_deion_O2p_Om__O2X_O3P_data;
% num_eV=1;           % num of levels up to wich eV processes are included
% R(1:num_eV)=R(1:num_eV)+R_eV_data(1:num_eV);
R(num_O)=-2*(sum(R_diss_data+R_eD)+R_rec_O2X_data+R_rec_O2a_data...
    +R_rec_O2b_data+R_rec_3O3P__O2X_O3P_data)+...
    R_mix_1DX_3Pb_data+R_deex_1D3PbyO2_data-R_ex_O3P_O1D_data+...  % O(3P)
    R_d_O2X_3P1D_data+R_d_O2a_3P1D_data...
    -R_rec_O3P_O2X_O3P__O3X_O3P_data-R_rec_O3P_O2X_O2__O3X_O2_data...
    +R_d_O2b_O3X__2O2X_O3P_data+R_d_O2a_O3X__2O2X_O3P_data...
    -R_exch_O3P_O3X__2O2X_data-R_exch_O3P_O3X__O2a_O2X_data...
    -R_exch_O3P_O3X__O2b_O2X_data+R_d_O2b_O3P_O1D_data...
    +R_deex_O1D_O3P__2O3P_data+R_mix_1DX_3Pa_data+R_d_O3X__O3P_O2X_data...
    -R_rec_O3P_O2X_O3X__2O3X_data+2*R_exch_O1D_O3X__O2X_2O3P_data...
    +2*R_diss_e_O2p__2O3P_data+R_diss_e_O2p__O3P_O1D_data...
    +R_exch_Op_O2X__O2p_O3P_data+R_exch_Op_O2a__O2p_O3P_data...
    -R_ion_O3P_Op_data+R_mix_O2X__e_O3P_Op_data...
    +R_mix_O2a__e_O3P_Op_data-R_mix_Om_O3P__e_O2X_data...
    +R_mix_Om_O2b__e_O3P_O2X_data+2*R_deion_Op_Om__2O3P_data...
    +R_deion_O2p_Om__O2X_O3P_data+R_mix_e_O2X__Om_O3P_data...
    +R_mix_e_O2a__Om_O3P_data+R_deion_e_Om__2e_O3P_data...
    -R_rec_O3P_2O2X__O3exc_O2X_data+R_mix_O2a_O3exc__2O2X_O3P_data...
    -R_exch_O3P_O3exc__2O2X_data;
R(num_O1D)=-R_mix_1DX_3Pb_data-R_deex_1D3PbyO2_data+...            % O(1D)
    R_ex_O3P_O1D_data+R_d_O2X_3P1D_data+R_d_O2a_3P1D_data...
    +R_d_O2b_O3P_O1D_data-R_deex_O1D_O3P__2O3P_data...
    -R_mix_1DX_3Pa_data-R_exch_O1D_O3X__2O2X_data...
    -R_exch_O1D_O3X__O2X_2O3P_data+R_diss_e_O2p__O3P_O1D_data;
R(num_Op)=-R_exch_Op_O3X__O2p_O2X_data-R_exch_Op_O2X__O2p_O3P_data...
    -R_exch_Op_O2a__O2p_O3P_data+R_ion_O3P_Op_data...              % O+
    +R_mix_O2X__e_O3P_Op_data+R_mix_O2a__e_O3P_Op_data...
    -R_deion_Op_Om__2O3P_data;
R(num_Om)=-R_mix_Om_O2a__e_O3X_data-R_mix_Om_O3P__e_O2X_data...    % O-
    -R_mix_Om_O2X__e_O3X_data-R_mix_Om_O2b__e_O3P_O2X_data...
    -R_deion_Op_Om__2O3P_data-R_deion_O2p_Om__O2X_O3P_data...
    +R_mix_e_O2X__Om_O3P_data+R_mix_e_O2a__Om_O3P_data...
    -R_deion_e_Om__2e_O3P_data;
R(num_O3)=R_rec_O3P_O2X_O3P__O3X_O3P_data...                       % O3
    +R_rec_O3P_O2X_O2__O3X_O2_data-R_d_O2b_O3X__2O2X_O3P_data...
    -R_d_O2a_O3X__2O2X_O3P_data-R_exch_O3P_O3X__2O2X_data...
    -R_exch_O3P_O3X__O2a_O2X_data-R_exch_O3P_O3X__O2b_O2X_data...
    -R_d_O3X__O3P_O2X_data+R_rec_O3P_O2X_O3X__2O3X_data...
    -R_exch_O1D_O3X__2O2X_data-R_exch_O1D_O3X__O2X_2O3P_data...
    -R_exch_Op_O3X__O2p_O2X_data+R_mix_Om_O2a__e_O3X_data...
    +R_mix_Om_O2X__e_O3X_data+R_deex_O3exc_O3P__O3X_O3P_data...
    +R_deex_O3exc_O2X__O3X_O2X_data;
R(num_O3exc)=R_rec_O3P_2O2X__O3exc_O2X_data...                     % O3exc
    -R_deex_O3exc_O3P__O3X_O3P_data-R_deex_O3exc_O2X__O3X_O2X_data...
    -R_mix_O2a_O3exc__2O2X_O3P_data-R_exch_O3P_O3exc__2O2X_data;
                                                                   
R=R*n0*tau;

% переделываем под уравнение для T
% cp=28.2+6456.2/788.3/sqrt(pi/2)*exp(-2*((T-1006.9)/788.3)^2);
% eV_to_J=1.60218e-19;
% O2a_diss=0;
% R_mm=0;
%     % R*dE
%     % O2(a) +O2(X)->2O(3P)+O2(X), O2(a) +O2(a)->2O(3P)+O2(a) 
%     %  and O2(a) +O(3P)->2O(3P)+O(3P), for rec 2O(3P)+O2(X)->O2(a) +O2(X),
%     %   2O(3P)+O2(a)->O2(a) +O2(a)
% % Q_H=(R_diss_data(num_O2a)-R_rec_O2a_data)*(O2.e_E(2)-2*O.form_e);
% Qin=(Q_el+0)*eV_to_J*ne_b*nAll_b*tau*n0/(k*T);
% R(num_T)=-Qin/(nAll_b*cp/N_a/k);

    out=R;
if mod(round(cputime), 100)==1
    if t==0
        disp_per=0;
    else
       disp_per=round(log10(log10(t+1)+1)/log10(log10(xspan(2)+1)+1)*100);
    end
    disp([num2str(t*tau) ' ' ...
        num2str(disp_per) '% ' num2str(T,'%4.12f')])
%     disp(['Rparts ' num2str(R(1)) ' ' num2str(R(num_T)) ' ' ...
%         num2str(R(num_O2a))])
end
% if mod(round(cputime), 10)==1
    pr_vec=[t*tau y(1:num_e)'*n0 y(num_T)*T0];
    fprintf(fileData, [repmat('%d ', 1, numel(pr_vec)), '\n'], pr_vec);
% end
% fclose('all');
% pause
% disp(t*tau)
end

function out=R_eV(n1, nO2a, nO2b, ne, nO, nO1D, nOm, nO3, EdN, ...
                                                            fileKrate, t, n0)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
    nAll=sum(n1)+nO2a+nO2b+nO+nO1D;
%     oxygen_saver(nAll*n0*k*T, T, ne*n0, n1*n0, ...
%         nO2a*n0, nO2b*n0, nO*n0, nO1D*n0, T_e, EdN);
%     cd LoKI-master\Code
% %     data4me=lokibcl_f('Oxygen\oxygen.in');
%     lokibcl_f('Oxygen\oxygen.in');
%     cd ..\..
    load LoKI-master\Code\Output\RNF\data4Maksim.mat
    data4me=data4Maksim;
    R(1:26)=0;
    data_length=length(data4me.rate);
    data_length_extra=length(data4me.rateExtra);
    
%     for ii=1:data_length
%         disp(ii)
%         disp(data4me.rate(ii).collDescription)
%     end
%     disp('EXTRA')
%     for ii=1:data_length_extra
%         disp(ii)
%         disp(data4me.rateExtra(ii).collDescription)
%     end

% if length(data4me.rate(8).value)~=2
%     save test.mat data4me
% end
% if convertCharsToStrings(data4me.rate(end).collDescription)~=...
%                                            "e+O(3P)<->e+O(1D),Excitation"
%     disp('WARNING e rates')
%     save test.mat data4me
% end
for i=1:data_length
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                           "e+O2(X)->e+2O(3P),Excitation"
        k_O2X_2O3P=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                        "e+O2(X)<->e+O2(a1Dg),Excitation"
        k_O2X_O2a=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                       "e+O2(X)<->e+O2(b1Sg+),Excitation"
        k_O2X_O2b=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                           "e+O(3P)<->e+O(1D),Excitation"
        k_O3P_O1D=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                      "e+O2(X)->e+O(3P)+O(1D),Excitation"
        k_O2X__O3P_O1D=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                   "e+O2(X)->e+e+O2(+,X),Ionization"
        k_O2X__e_O2p=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                   "e+O(3P)->e+e+O(+,gnd),Ionization"
        k_O3P__e_Op=data4me.rate(i).value;
    end
    if convertCharsToStrings(data4me.rate(i).collDescription)==...
                                   "e+O2(X)->O(-,gnd)+O(3P),Attachment"
        k_e_O2X__Om_O3P=data4me.rate(i).value;
    end
end
for i=1:data_length_extra
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                    "e+O2(a1Dg)<->e+O2(b1Sg+),Excitation"
        k_O2a_O2b=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O2(a1Dg)->e+O(3P)+O(1D),Excitation"
        k_O2a__O3P_O1D=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O2(a1Dg)->e+2O(3P),Excitation"
        k_O2a__2O3P=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O2(b1Sg+)->e+2O(3P),Excitation"
        k_O2b__2O3P=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                               "e+O2(b1Sg+)->e+O(3P)+O(1D),Excitation"
        k_O2b__O3P_O1D=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O3(X)->e+O(3P)+O2(X),Excitation"
        k_O3X__O3P_O2X=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O2(a1Dg)->e+e+O2(+,X),Ionization"
        k_O2a__e_O2p=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                               "e+O2(X)->e+e+O(3P)+O(+,gnd),Ionization"
        k_O2X__e_O3P_Op=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                            "e+O2(a1Dg)->e+e+O(3P)+O(+,gnd),Ionization"
        k_O2a__e_O3P_Op=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                                   "e+O2(a1Dg)->O(-,gnd)+O(3P),Attachment"
        k_e_O2a__Om_O3P=data4me.rateExtra(i).value;
    end
    if convertCharsToStrings(data4me.rateExtra(i).collDescription)==...
                            "e+O(-,gnd)->e+O(3P),Excitation"
        k_Om__e_O3P=data4me.rateExtra(i).value;
    end
end
% disp("krate")
% disp(k_Om__e_O3P)
    % e + O2(X) ->  e + 2O(3P)
R(6)=k_O2X_2O3P;
    % e + O2(X) <-> e + O2(a)
R(7)=ne*(sum(n1)*k_O2X_O2a(1)-nO2a*k_O2X_O2a(2));
% disp('krate')
% disp(k_O2X_O2a(1))
% disp(k_O2X_O2a(2))
    % e + O2(X) <-> e + O2(b)
R(8)=ne*(sum(n1)*k_O2X_O2b(1)-nO2b*k_O2X_O2b(2));
    % e + O2(a) <-> e + O2(b)
R(9)=ne*(nO2a*k_O2a_O2b(1)-nO2b*k_O2a_O2b(2));
    % e + O(3P) <-> e + O(1D)
R(10)=ne*(nO*k_O3P_O1D(1)-nO1D*k_O3P_O1D(2));
    % e + O2(X) -> e + O(3P) + O(1D)
R(11)=ne*(n1*k_O2X__O3P_O1D(1));
    % e + O2(a1Dg) -> e + O(3P) + O(1D)
R(12)=ne*(nO2a*k_O2a__O3P_O1D(1));%-...
            %nO*nO1D*n0*data4me.rateExtra(2).value(2));
    % e + O2(a) ->  e + 2O(3P)
R(13)=k_O2a__2O3P;
    % e + O2(b) ->  e + 2O(3P)
R(14)=k_O2b__2O3P;
    % e + O2(b) ->  e + O(3P) + O(1D)
R(15)=ne*nO2b*k_O2b__O3P_O1D;
    % e + O3(X) ->  e + O(3P) + O2(X)
R(16)=ne*nO3*k_O3X__O3P_O2X;
    % e + O2(X) -> 2e + O2(+,X)
R(17)=ne*n1*k_O2X__e_O2p;
    % e + O2(a1Dg) -> 2e + O2(+,X)
R(18)=ne*nO2a*k_O2a__e_O2p;
    % e + O(3P) -> 2e + O(+,gnd)
R(19)=ne*nO*k_O3P__e_Op;
    % e + O2(X) -> 2e + O(3P) + O(+,gnd)
R(20)=ne*n1*k_O2X__e_O3P_Op;
    % e + O2(a1Dg) -> 2e + O(3P) + O(+,gnd)
R(21)=ne*nO2a*k_O2a__e_O3P_Op;
    % e + O2(X) -> O(-,gnd) + O(3P)
R(22)=ne*n1*k_e_O2X__Om_O3P;
    % e + O2(a1Dg) -> O(-,gnd) + O(3P)
R(23)=ne*nO2a*k_e_O2a__Om_O3P;
    % e + O(-,gnd) -> 2e + O(3P)
R(24)=ne*nOm*k_Om__e_O3P;

R(end-1)=data4me.Te;
R(end)=data4me.Qin;
% pr_vec=[t k_O2X_2O3P k_O2X_O2a(1) k_O2X_O2a(2) k_O2X_O2b(1) k_O2X_O2b(2)...
%     k_O2a_O2b(1) k_O2a_O2b(2) k_O3P_O1D(1) k_O3P_O1D(2) k_O2X__O3P_O1D...
%     k_O2a__O3P_O1D];
% fprintf(fileKrate, [repmat('%d ', 1, numel(pr_vec)), '\n'], pr_vec);
if mod(round(cputime), 10)==1
    pr_vec(1:1+length(data4me.rate)+length(data4me.rateExtra))=0;
    pr_vec(1)=t;
    for i=1:length(data4me.rate)
        pr_vec(1+i)=data4me.rate(i).value(1);
    end
    for i=1:length(data4me.rateExtra)
        pr_vec(1+length(data4me.rate)+i)=data4me.rateExtra(i).value(1);
    end
    fprintf(fileKrate, [repmat('%d ', 1, numel(pr_vec)), '\n'], pr_vec);
end

% dE=M.e_i(1, 2:5)-M.e_i(1, 1);           % считаем k_down по детальному
% k0i_down=k0i_up.*k_bf_VT(T, -dE)';    % балансу
% gg=k_bf_VT(T, -dE);
% disp(gg(1))

% R(1)=(n1(2:5)'*k0i_down-sum(n1(1)*k0i_up))*ne;
% R(2:5)=(n1(1)*k0i_up-n1(2:5).*k0i_down)'*ne;
% R(6)=k_diss;
% R(7)=k_el;
% R(8)=M.e_i(1,2)*(k0i_up(1)*n1(1)-k0i_down(1)*n1(2))*ne+...     1
%      M.e_i(1,3)*(k0i_up(2)*n1(1)-k0i_down(2)*n1(3))*ne+...     2
%      M.e_i(1,4)*(k0i_up(3)*n1(1)-k0i_down(3)*n1(4))*ne+...     3
%      M.e_i(1,5)*(k0i_up(4)*n1(1)-k0i_down(4)*n1(5))*ne;      % 4
out=R';
% fileID = fopen('LoKI-master\Code\Output\RNF\eedf.txt','r');
% fgetl(fileID);
% for i_chk=2:500
%     str=fgetl(fileID);
%     temp=sscanf(str,'%f');
%     k_minus=temp(2);
%     if k_minus<0
%        disp('ERROR EEDF<0')
%        disp(k_minus)
%        disp(i_chk)
%        break
%     end
%     if k_minus==0
%        break
%     end
% end
% fclose(fileID);
end

function out = R_diss_e(T, ni_b, nO_b, ne_b, M, A1, A2, kd_eq, n0)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
h = 6.626070041e-34;    % постоянная Планка (Дж/с)
c = 299792458;
e_i=M.e_i(1,:);
% equil. vibr. partition function
ZvT = sum(exp(-e_i/(T*k)));
% parameter of TM model
U = Inf;
ZvU = sum(exp(e_i/(U*k)));
% non-equilibrium factor
Z = ZvT / ZvU * exp(e_i/k*(1/T + 1/U));
% dis. rates
kd = (kd_eq * Z)'; % m^3/sec

sigma = 2;      % 2 для гомоядерных, 1 для гетероядерных
Be=M.Be(1);%193.128087;
Theta_r = Be*h*c/k;
Z_rot = T./(sigma.*Theta_r);
Kdr2=(M.mass/A1.mass/A2.mass)^(3/2)*h^3*(2*pi*k*T)^(-3/2)*Z_rot...
    *exp(-(M.e_i(1,:)'+(M.form_e-A1.form_e-A2.form_e))/k/T);
kr=0;% kd .* Kdr2 * n0;     % всё говорит о том, что обратные  
                                % процессы не идут
RR = ne_b * nO_b*nO_b*kr;
RD = -ne_b*ni_b.*kd;
RD(end+1)=(M.diss_e-M.e_i(1,:)-M.e_0)*RD;
out = RD;
end

function out=R_VV(M1, M2, T, n1, n2) %VV for O2-O2 collision
k = 1.3807e-23; % постоянная Больцмана (Дж/К)

core1=2:M1.num_vibr_levels(1)-1; %c 1 по n-1 уровень
% core2=2:M2.num_vibr_levels-1; %c 1 по n-1 уровень

R=zeros(M1.num_vibr_levels(1), 1);
k_down=kvv_fho(T, M1);
dE1=M1.e_i(1,1:end-1)-M1.e_i(1,2:end);
dE2=M2.e_i(1,2:end)-M2.e_i(1,1:end-1);
kbf=exp((dE1'+dE2)/k/T);

k_up=k_down.*kbf;

R(1)=n1(2)*(k_down(1,:)*n2(1:end-1))-n1(1)*(k_up(1,:)*n2(2:end));
R(end)=n1(end-1)*(k_up(end,:)*n2(2:end))-...
    n1(end)*(k_down(end,:)*n2(1:end-1));
R(core1)=n1(core1+1).*(k_down(core1,:)*n2(1:end-1))...
    +n1(core1-1).*(k_up(core1-1,:)*n2(2:end))...
    -n1(core1).*(k_up(core1,:)*n2(2:end)+k_down(core1-1,:)*n2(1:end-1));

out=R;
end

function out=R_VT(M, T, n1, na)

core=2:M.num_vibr_levels(1)-1; %c 1 по n-1 уровень

R=zeros(M.num_vibr_levels(1), 1);
% k_down_O2=k_VT(M, T, 1:M.num_vibr_levels-1, -1, M)';
kvt_temp=kvt_fho(T, M);
k_down_O2=kvt_temp(1,:)';
k_down_O=kvt_temp(2,:)';
% k_up_O2=k_VT(M, T, 0:M.num_vibr_levels-2, 1, M)';
i=0:M.num_vibr_levels(1)-2;
delta_i=1;
dE=M.e_i(1,i+1)-M.e_i(1,i+1+delta_i);
k_up_O2=k_down_O2.*k_bf_VT(T, -dE)';
k_up_O=k_down_O.*k_bf_VT(T, -dE)';
% gg=k_bf_VT(T, -dE);
% disp(['VT ' num2str(gg(1)) ' ' num2str(1/gg(1))])

n2=sum(n1);

% заполняем серединку
R(core)=n2*(k_down_O2(core).*n1(core+1)+k_up_O2(core-1).*n1(core-1)-...
        k_down_O2(core-1).*n1(core)-k_up_O2(core).*n1(core))...
        +na*(k_down_O(core).*n1(core+1)+k_up_O(core-1).*n1(core-1)-...
        k_down_O(core-1).*n1(core)-k_up_O(core).*n1(core));
    
% заполняем для первого и последнего уровней
R(1)=n2*(k_down_O2(1)*n1(2)-k_up_O2(1)*n1(1))...
    +na*(k_down_O(1)*n1(2)-k_up_O(1)*n1(1));
R(M.num_vibr_levels(1))=...
    n2*( k_up_O2(M.num_vibr_levels(1)-1)*n1(M.num_vibr_levels(1)-1)-...
    k_down_O2(M.num_vibr_levels(1)-1)*n1(M.num_vibr_levels(1)) )...
    +na*( k_up_O(M.num_vibr_levels(1)-1)*n1(M.num_vibr_levels(1)-1)-...
    k_down_O(M.num_vibr_levels(1)-1)*n1(M.num_vibr_levels(1)) );
out=R;
end

function out=k_VT(M, T, i, delta_i, N)
% это через интегралы, надо бы доделать и внедрить
HBAR = 1.054571800e-34;
alpha_FHO =4e10;
E_FHO =200;
SVT_FHO=0.44;
    
    dE=M.e_i(1,i+1)-M.e_i(1,i+1+delta_i);
    coll_mass=M.mass*N.mass/(M.mass+N.mass);
    coll_diameter=0.5*(M.diameter+N.diameter);
	if delta_i < 0
		omega = - dE / (HBAR * delta_i);
		res = k_VT_FHO_RS(T, coll_mass, coll_diameter, M.red_osc_mass,...
            dE, i, delta_i, omega, M.mA, alpha_FHO, E_FHO, SVT_FHO);
		out=0.5 * (res + k_VT_FHO_RS(T, coll_mass, coll_diameter, ...
            M.red_osc_mass, dE, i, delta_i, omega, M.mB, ...
            alpha_FHO, E_FHO, SVT_FHO));
	else 
		omega = -dE / (HBAR* delta_i);
		res = k_VT_FHO_RS(T, coll_mass, coll_diameter, M.red_osc_mass, ...
            -dE, i + delta_i, -delta_i, omega, M.mA, alpha_FHO, E_FHO, ...
            SVT_FHO);
		res = 0.5 * (res + k_VT_FHO_RS(T, coll_mass, coll_diameter, ...
            M.red_osc_mass, -dE, i + delta_i, -delta_i, omega, M.mB, ...
            alpha_FHO, E_FHO, SVT_FHO));
		out=res .* k_bf_VT(T, -dE);
    end
end

function out=k_VT_FHO_RS(T, coll_mass, diameter, osc_mass, ...
    Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO)
	out=8 * integral_VT_FHO_RS(T, 0, coll_mass, diameter, osc_mass, ...
        Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO);
end

function out=k_bf_VT(T, dE) 
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
	out=exp(-dE / (k * T));
end

function out=integral_VT_FHO_RS(T, degree, coll_mass, diameter, ...
    osc_mass, Delta_E_vibr, j, delta_i, omega, ram,...
	alpha_FHO, E_FHO, svt_FHO)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
	conversion = sqrt(2 * k * T / coll_mass);
    integrant=@(g) g.^(2*degree+3).*...
        crosssection_VT_FHO_RS(conversion.*g,coll_mass, diameter, ...
        osc_mass, Delta_E_vibr, j, delta_i, omega, ram, alpha_FHO, ...
        E_FHO, svt_FHO) .* exp(-g .* g);
	out = sqrt(k * T / (2 * pi * coll_mass)) * integral(integrant, 0,...
        inf, 'ArrayValued', true);
end

function out=crosssection_VT_FHO_RS(rel_vel, coll_mass, diameter,...
    osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO,...
    svt_FHO)
	if delta_i < 0
		out=probability_VT_FHO(rel_vel, coll_mass, osc_mass, ...
            Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, ...
            svt_FHO) * crosssection_elastic_RS(diameter);
    else
		out = 0;
    end
end

function out=probability_VT_FHO(rel_vel, coll_mass, osc_mass, ...
    Delta_E_vibr, i, delta_i,	omega, ram, alpha_FHO, E_FHO, svt_FHO)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
h = 6.626070041e-34;    % постоянная Планка (Дж/с)
    avg_vel = vel_avg_vt(rel_vel, coll_mass, Delta_E_vibr);
    if delta_i > 0
        s = delta_i;
        ns = factorial(i + delta_i) ./ factorial(i);
    else
        s = -delta_i;
        ns = factorial(i) ./ factorial(i + delta_i);
    end
    phi = (2. / pi) * atan(sqrt((2 * k* E_FHO) ./ (coll_mass * ...
        avg_vel .* avg_vel)));
    eps = cosh((1 + phi) .* pi .* omega ./ (alpha_FHO * avg_vel));
    eps = eps .* (pi .* ram .* coll_mass ./ (sinh(2 * pi .* omega ./...
        (alpha_FHO .* avg_vel)) * alpha_FHO));
    eps = eps .* eps;
    eps = eps .*(8 * svt_FHO * omega / (osc_mass * h));
    out = ns .* eps.^s .* exp(-2.0 .* ns.^ (1.0 ./ s).* eps ./...
        (s + 1.0)) ./ (factorial(s) .* factorial(s));
    out(avg_vel <= 0.0)=0;
end

function out=crosssection_elastic_RS(diameter)
	out = pi * diameter * diameter;
end

function out=vel_avg_vt(rel_vel, coll_mass, Delta_E_vibr)
    rel_vel_after_sq = Delta_E_vibr * (2.0 / coll_mass) + ...
        rel_vel .* rel_vel;
    out = 0.5 * (rel_vel + sqrt(rel_vel_after_sq));
    out(out<0)=-1;
end

function out = R_diss(T, ni_b, nO_b, n0, M, A1, A2)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
h = 6.626070041e-34;    % постоянная Планка (Дж/с)
c = 299792458;

CA=[0.0330    0.1659]*1e-7;
nA=-1.5;
D=59379;
e_i=M.e_i(1,:);
% equil. coef-s
kd_eq = CA.*T.^nA*exp(-D/T);% m^3/sec
% equil. vibr. partition function
ZvT = sum(exp(-e_i/(T*k)));
% parameter of TM model
U = Inf;

ZvU = sum(exp(e_i/(U*k)));

% non-equilibrium factor
Z = ZvT / ZvU * exp(e_i/k*(1/T + 1/U));

% dis. rates
kd = (kd_eq' * Z)'; % m^3/sec
kd(2:end)=0;    % turn on back later for higher levels

sigma = 2;      % 2 для гомоядерных, 1 для гетероядерных
Be=M.Be(1);%193.128087;
Theta_r = Be*h*c/k;
Z_rot = T./(sigma.*Theta_r);
Kdr2=(M.mass/A1.mass/A2.mass)^(3/2)*h^3*(2*pi*k*T)^(-3/2)*Z_rot...
    *exp(-(M.e_i(1,:)'+(M.form_e-A1.form_e-A2.form_e))/k/T);
kr= kd .* Kdr2 * n0;
kd_O2=kd(:,1);
kd_O=kd(:,2);
kr_O2=kr(:,1);
kr_O=kr(:,2);

nm_b=sum(ni_b);

RD = nm_b * (nO_b*nO_b*kr_O2-ni_b.*kd_O2) + ...
        nO_b * (nO_b*nO_b*kr_O-ni_b.*kd_O);
out = RD;
end