function out=Rpart_Tiago_STS(t, y, Prcl, Coll, tau, xspan, ...
                        T0, n0, is_e, EdN, fileData, fileKrate, config)
% правые части для системы ДУ для нульмерной релаксации в кислороде с 
% электронным возбуждением. Хотим повторить результаты Тиаго, поэтому
% делаем закос под его код.
% 26.11.2020
k = 1.3807e-23;  % постоянная Больцмана (Дж/К)
N_a = 6.0221e23; % постоянная Авогадро (1/моль)

O2=Prcl.O2; O=Prcl.O;
num_O2=1;%O2.num_vibr_levels(1);
num_O2a=num_O2+O2.num_vibr_levels(1);
num_O2b=num_O2a+O2.num_vibr_levels(2);
num_O2p=num_O2b+O2.num_vibr_levels(3);
num_O=num_O2p+1;
num_O1D=num_O+1;
num_Op=num_O1D+1;
num_Om=num_Op+1;
num_O3=num_Om+1;
num_O3exc=num_O3+1;
num_e=num_O3exc+1;
num_T=num_e+1;
R=zeros(num_T,1);
nO2_b=y(num_O2:num_O2-1+O2.num_vibr_levels(1));
nm=sum(y(1:num_O2b-1+O2.num_vibr_levels(3)));%sum(nO2_b);
nO_b=y(num_O);
ne_b=y(num_e);
T=y(num_T)*T0;
nAll_b=sum(y(1:num_O3exc));
if nAll_b~=nm+y(num_O2p)+nO_b+y(num_O1D)+y(num_Op)+y(num_Om)+y(num_O3)+...
        y(num_O3exc)
    disp('ALERT')
    disp(t*tau)
    disp(T)
    disp(ne_b*n0)
    disp(y(num_O2p))
    pause
end

R_VT_data=zeros(num_O2-1+O2.num_vibr_levels(1), 1);
R_VV_data=zeros(num_O2-1+O2.num_vibr_levels(1), 1);
R_diss_data=zeros(num_O2b-1+O2.num_vibr_levels(3),1);
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
R_O2Xi_O2a_data=0;  R_O2Xi_O2b_data=0;  
R_O2Xi_2O3P_data=zeros(O2.num_vibr_levels(1),1);
R_O2Xi__O3P_O1D_data=zeros(O2.num_vibr_levels(1),1);
R_O2Xi__O3P_Om_data=zeros(O2.num_vibr_levels(1),1);
R_O2Xi_e_O2p_data=zeros(O2.num_vibr_levels(1),1);
R_eV_cs_data=0;

Qin=0;
    % vibrational kinetics
% надо бы тут VT от разных других включить, а не только O2(X) и O(3P)
switch config.V
    case "Paper"
        R_VT_data=Rvt_O2O2_Ann(nO2_b, T, O2)...
                        +Rvt_O2O_Esposito(nO2_b, nO_b+y(num_O1D), T, O2);
        R_VV_data=Rvv_O2O2_Ann(nO2_b, T, O2);
    case "FHO"
        R_VT_data=R_VT(O2, T, nO2_b, nO_b+y(num_O1D), config.V);
        R_VV_data=R_VV(O2, O2, T, nO2_b, nO2_b);
    case "FHO10"
        R_VT_data=R_VT(O2, T, nO2_b, nO_b+y(num_O1D), config.V);
        R_VV_data=R_VV(O2, O2, T, nO2_b, nO2_b);
end
%     R_VT_data=R_VT(O2, T, nO2_b, nO_b+y(num_O1D));
%     R_VV_data=R_VV(O2, O2, T, nO2_b, nO2_b);
%     R_VT_data=Rvt_O2O2_Ann(nO2_b, T, O2)...
%                         +Rvt_O2O_Esposito(nO2_b, nO_b+y(num_O1D), T, O2);
%     R_VV_data=Rvv_O2O2_Ann(nO2_b, T, O2);

%     if nO_b>0
%     gg=Rvv_O2O2_Ann(nO2_b, T, O2);
% %     disp(R_VV_data)
% %     disp(gg)
%     gg=[R_VV_data gg];
%     disp(gg)
%         pause
%     end

%     disp(R_VT_data)
%     if nO_b>0
%     gg=Rvt_O2O_Esposito(nO2_b, nO_b, T, O2);
%     gg=[R_VT_data gg];
%     disp(gg)
%         pause
%     end

%     if nO_b>0
%     gg=Rvt_O2O2_Ann(nO2_b, T, O2);
%     gg=[R_VT_data gg];
%     disp(gg)
%         pause
%     end

    
    % O2(X,a,b)+O2(X,a,b)->2O(3P)+O2(same)
% ! дисскоциация от тяжёлых частиц пока отключена, потому что в задачах с
%   разрядами её не рассматривают. Её нужно протестить будет для
%   корректности с колебательной кинетикой
% ! Сейчас подключена модель МТ, потому что Тиаго использует диссоциацию.
%   Но вообще нужно с неё уйти на модель португальцев или на Алиата
% R_diss_data=R_diss_data + R_diss_Aliat(... 
%     Prcl.O2, y(1:num_O2b), nO_b, nO_b, nm, Coll.O2_O2, T,n0);
U=Inf;
ZvT = sum(exp(-O2.e_i(1,:)/(T*k)));
ZvU = sum(exp(O2.e_i(1,:)/(U*k)));
Z = ZvT / ZvU * exp(O2.e_i(1,:)/k*(1/T + 1/U));
    kd_eqO2O2 = Coll.O2_O2.ArrA(1) * T^Coll.O2_O2.ArrN(1)*...
                                                exp(-O2.diss_e(1)/(k*T)); 
    kdO2O2 = kd_eqO2O2' * Z'; % m^3/sec
    R_O2O2=kdO2O2.*nO2_b*sum(nO2_b);
R_diss_data(1:42)=R_diss_data(1:42) - R_O2O2;
    % O2(X,a,b)+O(3P)->2O(3P)+O(3P)
% R_diss_data=R_diss_data + R_diss_Aliat(Prcl.O2, ...
%     y(1:num_O2b), nO_b, nO_b, sum(y(num_O:num_Om)), Coll.O2_O, T, n0);
    kd_eqO2O = Coll.O2_O.ArrA(1) * T^Coll.O2_O.ArrN(1)*...
                                                exp(-O2.diss_e(1)/(k*T)); 
    kdO2O = kd_eqO2O' * Z'; % m^3/sec
    R_O2O=kdO2O.*nO2_b*(y(num_O)+y(num_O1D));
R_diss_data(1:42)=R_diss_data(1:42) - R_O2O;
    % Вводим Diss для O2a и O2b, потому что чтобы хотя бы какая-то была
R_diss_data(num_O2a)=R_diss_data(num_O2a)...
                            -kd_eqO2O2*y(num_O2a)*sum(nO2_b)...
                            -kd_eqO2O*y(num_O2a)*(y(num_O)+y(num_O1D));
R_diss_data(num_O2b)=R_diss_data(num_O2b)...
                            -kd_eqO2O2*y(num_O2b)*sum(nO2_b)...
                            -kd_eqO2O*y(num_O2b)*(y(num_O)+y(num_O1D));

% if is_e~=0
%     if t*tau>1e-3
%         run_LoKI(T, nO2_b, y(num_O2a), y(num_O2b), ...
%             ne_b, nO_b, y(num_O1D), n0, T_e, EdN, t*tau);
%     end
% end
if is_e~=0
    switch config.LoKI
        case 0
            run_LoKI_STS(T, nO2_b, y(num_O2a), y(num_O2b), ...
                                ne_b, nO_b, y(num_O1D), n0, EdN, t*tau);
        case -1
        otherwise
            if mod(round(cputime), 8)==1
                run_LoKI_STS(T, nO2_b, y(num_O2a), y(num_O2b), ...
                                ne_b, nO_b, y(num_O1D), n0, EdN, t*tau);
            end
    end
end

if is_e~=0
        % for MT model 
    U=Inf;
    ZvT = sum(exp(-O2.e_i(1,:)/(T*k)));
    ZvU = sum(exp(O2.e_i(1,:)/(U*k)));
    Z = ZvT / ZvU * exp(O2.e_i(1,:)/k*(1/T + 1/U));
%     kd_eqO2O2 = Coll.O2_O2.ArrA(1) * T^Coll.O2_O2.ArrN(1)*...
%                                                     exp(-O2.diss_e(1)/(k*T)); 
%     kd = kd_eqO2O2' * Z'; % m^3/sec
%     disp(num2str([(0:41)' kd]))
%     pause
        % call function
    [R_eV_data, R_O2Xi_O2a_data, R_O2Xi_O2b_data, R_O2Xi_2O3P_data,...
        R_O2Xi__O3P_O1D_data, R_O2Xi__O3P_Om_data, R_O2Xi_e_O2p_data,...
        R_eV_cs_data]=...
        R_eV(nO2_b, y(num_O2a), y(num_O2b), ne_b, nO_b, y(num_O1D), ...
        y(num_Om), y(num_O3), fileKrate, t*tau, Z');
end
    % посчитать kd по МТ и сравнить с Тиаго
keq=1;
    % recombination 2O(3P)+O2(X,0)->O2(X,0)+O2(X,0)%можно nm вместо nO2_b?
R_rec_O2X_data=...
           0.5*3.81e-30/T*exp(-170/T)/1e12*y(num_O)*y(num_O)*y(num_O2)*n0;
    % recombination 2O(3P)+O2(X,0)->O2(a)+O2(X,0) % можно nm вместо nO2_b?
R_rec_O2a_data=0.33*3.81e-30/T*exp(-170/T)/1e12*nO_b*nO_b*y(num_O2)*n0;
    % recombination 2O(3P)+O2(X,0)->O2(b)+O2(X,0) % можно nm вместо nO2_b?
R_rec_O2b_data=0.17*3.81e-30/T*exp(-170/T)/1e12*nO_b*nO_b*y(num_O2)*n0;
    % recombination 2O(3P)+O(3P)->O2(X,0)+O(3P)
R_rec_3O3P__O2X_O3P_data=...
                    3.6e-32*(1/T)^0.63/1e12*y(num_O)*y(num_O)*y(num_O)*n0;
    % k_eD
k_eD_eq=R_eV_data(6);
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
if is_e==1
        % diss e + O2(+,X) -> 2O(3P)            % not affect on ne
    R_diss_e_O2p__2O3P_data=2e-7*(300/Te)/1e6*ne_b*y(num_O2p);
        % diss e + O2(+,X) -> O(3P) + O(1D)     % not affect on ne
    R_diss_e_O2p__O3P_O1D_data=1.95e-7*(300/Te)^0.7/1e6*ne_b*y(num_O2p);
end
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
        % Other electron processes. Uncomment to exclude.
    % e + O2(X, 1-6) -> e + O2(a)
% R_O2Xi_O2a_data=0;
    % e + O2(X, 1-8) -> e + O2(b)
% R_O2Xi_O2b_data=0;
    % e + O2(X, 1-41) -> e + 2O(3P)
% R_O2Xi_2O3P_data=zeros(O2.num_vibr_levels(1),1);
    % e + O2(X, 1-41) -> e + O(3P) + O(1D)
% R_O2Xi__O3P_O1D_data=0;
    % e + O2(X, 1-41) -> O(3P) + O-
% R_O2Xi__O3P_Om_data=zeros(O2.num_vibr_levels(1),1);
    % e + O2(X, 1-41) -> 2e + O2+
% R_O2Xi_e_O2p_data=zeros(O2.num_vibr_levels(1),1);

    % deexcitation O2(b)+O(3P)->O2(X,0)+O(3P)
R_deex_O2b_O2_byO3p =4e-14/1e6*y(num_O2b)*nO_b;
    % deexcitation O2(b)+O(3P)->O2(a)+O(3P)
R_deex_O2b_O2a_byO3p=4e-14/1e6*y(num_O2b)*nO_b;
    % deexcitation O(1D)+O2(X,0)->O(3P)+O2(X,0)
R_deex_1D3PbyO2_data=7e-12*exp(67/T)/1e6*y(num_O1D)*y(num_O2);
    % deexcitation O2(a) + O(3P) -> O2(X,0) + O(3P)
R_deex_O2a_O3P_O2X_O3P_data=7e-17/1e6*y(num_O2a)*y(num_O);
    % deexcitation O2(a)+O2(X,0)->O2(X,0)+O2(X,0) % можно nm вместо nO2_b?
R_deex_O2a_O2__O2X_O2_data=2.2e-18*(T/300)^0.8/1e6*y(num_O2a)*y(num_O2);
    % deexcitation O(1D) + O(3P) -> O(3P) + O(3P)
R_deex_O1D_O3P__2O3P_data=8e-12/1e6*y(num_O1D)*y(num_O);

    % exchange O(1D) + O2(X,0) -> O(3P) + O2(b1Sg+)
R_mix_1DX_3Pb_data=2.56e-11*exp(67/T)/1e6*y(num_O1D)*y(num_O2);
    % exchange O(1D) + O2(X,0) -> O(3P) + O2(a1Dg)
R_mix_1DX_3Pa_data=1e-12/1e6*y(num_O1D)*y(num_O2);
    % exchange 2O2(a1Dg) -> O2(b1Sg+) + O2(X,0)
R_mix_2O2a_O2XO2b_data=...
               1.81e-18*exp(700/T)*(T/300)^3.8/1e6*y(num_O2a)*y(num_O2a);

    % O3
    % O(3P) + O2(X,0) + O(3P) -> O3(X) + O(3P)
R_rec_O3P_O2X_O3P__O3X_O3P_data=...
                           2.1e-34*exp(345/T)/1e12*y(num_O2)*nO_b*nO_b*n0;
    % O(3P) + O2(X,0) + O2(X,0) -> O3(X) + O2(X,0) %можно nm вместо nO2_b?
R_rec_O3P_O2X_O2__O3X_O2_data=...
                0.33*6.4e-35*exp(663/T)/1e12*nO_b*y(num_O2)*y(num_O2)*n0;
    % O2(b1Sg+) + O3(X) -> 2O2(X,0) + O(3P)
R_d_O2b_O3X__2O2X_O3P_data=1.5e-11/1e6*y(num_O2b)*y(num_O3);
    % O2(a1Dg) + O3(X) -> 2O2(X,0) + O(3P)
R_d_O2a_O3X__2O2X_O3P_data=...
                        5.2e-11*exp(-2840/T)/1e6*y(num_O2a)*y(num_O3);
    % O(3P) + O3(X) -> 2O2(X,0)
R_exch_O3P_O3X__2O2X_data=...
                        0.5*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
    % O(3P) + O3(X) -> O2(a1Dg) + O2(X,0)
R_exch_O3P_O3X__O2a_O2X_data=...
                        0.33*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
    % O(3P) + O3(X) -> O2(b1Sg+) + O2(X,0)
R_exch_O3P_O3X__O2b_O2X_data=...
                        0.17*1.8e-11*exp(-2300/T)/1e6*y(num_O)*y(num_O3);
	% O(3P) + O2(X,0) + O3(X) -> O3(X) + O3(X) (low)
R_rec_O3P_O2X_O3X__2O3X_data=...
                1.66e-34*exp(T/300)/1e12*y(num_O)*y(num_O2)*y(num_O3)*n0;
    % O(1D) + O3(X) -> 2O2(X,0) (low)
R_exch_O1D_O3X__2O2X_data=1.2e-10/1e6*y(num_O1D)*y(num_O3);
    % O(1D) + O3(X) -> O2(X,0) + 2O(3P) (low)
R_exch_O1D_O3X__O2X_2O3P_data=1.2e-10/1e6*y(num_O1D)*y(num_O3);

    % O+
    % exch O(+,gnd) + O3(X) -> O2(+,X) + O2(X,0)      % not affect on ne
R_exch_Op_O3X__O2p_O2X_data=1e-10/1e6*y(num_Op)*y(num_O3);
    % exch O(+,gnd) + O2(X) -> O2(+,X) + O(3P)      % not affect on ne
R_exch_Op_O2X__O2p_O3P_data=2e-11*(300/T)^0.5/1e6*y(num_Op)*y(num_O2);
    % exch O(+,gnd) + O2(a1Dg) -> O2(+,X) + O(3P)   % not affect on ne
R_exch_Op_O2a__O2p_O3P_data=2e-11*(300/T)^0.5/1e6*y(num_Op)*y(num_O2a);

    % O-
    % O(-,gnd) + O2(a1Dg) -> e + O3(X)
R_mix_Om_O2a__e_O3X_data=0.75*1.9e-10/1e6*y(num_Om)*y(num_O2a);
    % O(-,gnd) + O(3P) -> e + O2(X,0)  тут krate Тиаго, а не из статьи
R_mix_Om_O3P__e_O2X_data=1.9e-10/1e6*y(num_Om)*y(num_O);
    % O(-,gnd) + O2(X,0) -> e + O3(X)
R_mix_Om_O2X__e_O3X_data=1e-12/1e6*y(num_Om)*y(num_O2);
    % O(-,gnd) + O2(b1Sg+) -> e + O(3P) + O2(X,0)
R_mix_Om_O2b__e_O3P_O2X_data=6.9e-10/1e6*y(num_Om)*y(num_O2b);
    % O(+,gnd) + O(-,gnd) -> 2O(3P)
R_deion_Op_Om__2O3P_data=2.8e-7/1e6*y(num_Op)*y(num_Om);
    % O2(+,X) + O(-,gnd) -> O2(X,0) + O(3P)
R_deion_O2p_Om__O2X_O3P_data=9.6e-8*(300/T)^0.5/1e6*y(num_O2p)*y(num_Om);

    % O3(exc)
    % O(3P) + O2(X,0) + O2(X,0) -> O3(exc) + O2(X,0)
R_rec_O3P_2O2X__O3exc_O2X_data=...
            0.67*6.4e-35*exp(663/T)/1e12*y(num_O)*y(num_O2)*y(num_O2)*n0;
    % O3(exc) + O(3P) -> O3(X) + O(3P)
R_deex_O3exc_O3P__O3X_O3P_data=2e-13/1e6*y(num_O3exc)*y(num_O);
    % O3(exc) + O2(X,0) -> O3(X) + O2(X,0)
R_deex_O3exc_O2X__O3X_O2X_data=3e-15/1e6*y(num_O3exc)*y(num_O2);
    % O2(a1Dg) + O3(exc) -> 2O2(X) + O(3P)
R_mix_O2a_O3exc__2O2X_O3P_data=...
                26e-11*exp(-(2840-1553)/T)/1e6*y(num_O2a)*y(num_O3exc);
    % O(3P) + O3(exc) -> 2O2(X,0)
R_exch_O3P_O3exc__2O2X_data=...
                    8e-12*exp(-(2060-1553)/T)/1e6*y(num_O)*y(num_O3exc);
                 
	% O2(X,a,b) diss by Aliat model, не забыть перепроверить
R_eD=zeros(num_O2b-1+O2.num_vibr_levels(3),1);
% R_eD=R_diss_e_Aliat(Prcl.O2, ...
%     y(1:num_O2b-1+O2.num_vibr_levels(3)), nO_b, nO_b, ne_b, ...
%     [k_eD_eq k_eD_O2a_2O3P k_eD_O2b_2O3P], T, n0);
    % пока Алиата не использую. Мб позже.
    % e + O2(?) -> e + O(3P) + O(3P) by Phelps and Tashiro
R_eD(num_O2)=-k_eD_eq*ne_b*y(num_O2);
R_eD(num_O2a)=-k_eD_O2a_2O3P*ne_b*y(num_O2a);
R_eD(num_O2b)=-k_eD_O2b_2O3P*ne_b*y(num_O2b);
    % e + O2(X,v=1-41) -> e + 2O(3P) / e + O(3P) + O(1D)
R_eD(num_O2+1:num_O2-1+O2.num_vibr_levels(1))=...
                   -R_O2Xi_2O3P_data(2:end)-R_O2Xi__O3P_O1D_data(2:end);
        
%               ******собираем все R******
    % O2(X, i), VT, VV
R(num_O2:num_O2-1+O2.num_vibr_levels(1))=R_VT_data+R_VV_data;
    % O2(all), diss, eD
R(num_O2:num_O2b-1+O2.num_vibr_levels(3))=...
            R(num_O2:num_O2b-1+O2.num_vibr_levels(3))+R_diss_data+R_eD;
    % O2 eV processes
% num_eV=5;           % num of levels up to wich eV processes are included
% R(num_O2:num_O2-1+num_eV)=R(num_O2:num_O2-1+num_eV)+R_eV_data(1:num_eV);
R(num_O2:num_O2-1+O2.num_vibr_levels(1))=...
                    R(num_O2:num_O2-1+O2.num_vibr_levels(1))+R_eV_cs_data;
    % e + O2(X, 1-6) -> e + O2(a) / uncomment zero to exclude
R(num_O2:num_O2-1+6)=R(num_O2:num_O2-1+6)-R_O2Xi_O2a_data;
    % e + O2(X, 1-8) -> e + O2(b) / uncomment zero to exclude
R(num_O2:num_O2-1+8)=R(num_O2:num_O2-1+8)-R_O2Xi_O2b_data;
    % e + O2(X, 1-41) -> O(3P) + O- / uncomment zero to exclude 
R(num_O2+1:num_O2-1+O2.num_vibr_levels(1))=...          % not affect on ne
    R(num_O2+1:num_O2-1+O2.num_vibr_levels(1))-R_O2Xi__O3P_Om_data(2:end);
    % e + O2(X, 1-41) -> 2e + O2+ / uncomment zero to exclude 
R(num_O2+1:num_O2-1+O2.num_vibr_levels(1))=...          % not affect on ne
    R(num_O2+1:num_O2-1+O2.num_vibr_levels(1))-R_O2Xi_e_O2p_data(2:end);
    % aditional processes
R(num_O2)=R(num_O2)-R_ex_O2_O2a_data-R_ex_O2_O2b_data+...       % O2(X, 0)
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
    -R_mix_e_O2a__Om_O3P_data-R_mix_O2a_O3exc__2O2X_O3P_data...
    +sum(R_O2Xi_O2a_data);
R(num_O2b)=R(num_O2b)+R_rec_O2b_data+R_ex_O2_O2b_data...           % O2b
    -R_deex_O2b_O2_byO3p-R_deex_O2b_O2a_byO3p+R_ex_O2a_O2b_data...
    +R_mix_1DX_3Pb_data+R_mix_2O2a_O2XO2b_data...
    -R_d_O2b_O3X__2O2X_O3P_data+R_exch_O3P_O3X__O2b_O2X_data...
    -R_d_O2b_O3P_O1D_data-R_mix_Om_O2b__e_O3P_O2X_data...
    +sum(R_O2Xi_O2b_data);
R(num_O2p)=R_ion_O2X__e_O2p_data+R_ion_O2a__e_O2p_data...
    -R_diss_e_O2p__2O3P_data-R_diss_e_O2p__O3P_O1D_data...
    +R_exch_Op_O3X__O2p_O2X_data+R_exch_Op_O2X__O2p_O3P_data...
    +R_exch_Op_O2a__O2p_O3P_data-R_deion_O2p_Om__O2X_O3P_data...
    +sum(R_O2Xi_e_O2p_data);
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
    -R_exch_O3P_O3exc__2O2X_data+sum(R_O2Xi__O3P_Om_data);
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
    -R_deion_e_Om__2e_O3P_data+sum(R_O2Xi__O3P_Om_data);
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
    disp([num2str(t*tau, '%.15e') ' ' ...
        num2str(disp_per) '% ' num2str(T,'%4.12f')])
end
% if mod(round(cputime), 10)==1
%     pr_vec=[t*tau y(1:num_e)'*n0 y(num_T)*T0];
%     fprintf(fileData, [repmat('%d ', 1, numel(pr_vec)), '\n'], pr_vec);
% end
% fclose('all');
% pause
% disp(t*tau)
end

function [out, R_O2Xi_O2a, R_O2Xi_O2b, R_O2Xi_2O3P, R_O2Xi__O3P_O1D,...
                R_O2Xi__O3P_Om, R_O2Xi_e_O2p, R_eV_cs]=...
            R_eV(n1, nO2a, nO2b, ne, nO, nO1D, nOm, nO3, fileKrate, t, Z)
k = 1.3807e-23; % постоянная Больцмана (Дж/К)
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
%                                           "e+O(3P)<->e+O(1D),Excitation"
%     disp('WARNING e rates')
%     save test.mat data4me
% end

k_O2Xi_2O3P=zeros(1, 42);
k_eV_f=zeros(42);
k_eV_b=zeros(42);
check_arr=zeros(42);
for i=1:data_length
    switch convertCharsToStrings(data4me.rate(i).collDescription)
        case "e+O2(X,v=0)<->e+O2(X,v=1),Vibrational"
            k_eV(1, :)=data4me.rate(i).value;
        case "e+O2(X,v=0)<->e+O2(X,v=2),Vibrational"
            k_eV(2, :)=data4me.rate(i).value;
        case "e+O2(X,v=0)<->e+O2(X,v=3),Vibrational"
            k_eV(3, :)=data4me.rate(i).value;
        case "e+O2(X,v=0)<->e+O2(X,v=4),Vibrational"
            k_eV(4, :)=data4me.rate(i).value;
        case "e+O2(X)->e+2O(3P),Excitation"
            k_O2X_2O3P=data4me.rate(i).value;
        case "e+O2(X)<->e+O2(a1Dg),Excitation"
            k_O2X_O2a=data4me.rate(i).value;
        case "e+O2(X)<->e+O2(b1Sg+),Excitation"
            k_O2X_O2b=data4me.rate(i).value;
        case "e+O(3P)<->e+O(1D),Excitation"
            k_O3P_O1D=data4me.rate(i).value;
        case "e+O2(X)->e+O(3P)+O(1D),Excitation"
            k_O2X__O3P_O1D=data4me.rate(i).value;
        case "e+O2(X)->e+e+O2(+,X),Ionization"
            k_O2X__e_O2p=data4me.rate(i).value;
        case "e+O(3P)->e+e+O(+,gnd),Ionization"
            k_O3P__e_Op=data4me.rate(i).value;
        case "e+O2(X)->O(-,gnd)+O(3P),Attachment"
            k_e_O2X__Om_O3P=data4me.rate(i).value;
    end
end
for i=1:data_length_extra
 switch convertCharsToStrings(data4me.rateExtra(i).collDescription)
  case "e+O2(a1Dg)<->e+O2(b1Sg+),Excitation"
    k_O2a_O2b=data4me.rateExtra(i).value;
  case "e+O2(a1Dg)->e+O(3P)+O(1D),Excitation"
    k_O2a__O3P_O1D=data4me.rateExtra(i).value;
  case "e+O2(a1Dg)->e+2O(3P),Excitation"
    k_O2a__2O3P=data4me.rateExtra(i).value;
  case "e+O2(b1Sg+)->e+2O(3P),Excitation"
    k_O2b__2O3P=data4me.rateExtra(i).value;
  case "e+O2(b1Sg+)->e+O(3P)+O(1D),Excitation"
    k_O2b__O3P_O1D=data4me.rateExtra(i).value;
  case "e+O3(X)->e+O(3P)+O2(X,v=0),Excitation"
    k_O3X__O3P_O2X=data4me.rateExtra(i).value;
  case "e+O2(a1Dg)->e+e+O2(+,X),Ionization"
    k_O2a__e_O2p=data4me.rateExtra(i).value;
  case "e+O2(X,v=0)->e+e+O(3P)+O(+,gnd),Ionization"
    k_O2X__e_O3P_Op=data4me.rateExtra(i).value;
  case "e+O2(a1Dg)->e+e+O(3P)+O(+,gnd),Ionization"
    k_O2a__e_O3P_Op=data4me.rateExtra(i).value;
  case "e+O2(a1Dg)->O(-,gnd)+O(3P),Attachment"
    k_e_O2a__Om_O3P=data4me.rateExtra(i).value;
  case "e+O(-,gnd)->e+O(3P),Excitation"
    k_Om__e_O3P=data4me.rateExtra(i).value;
  case "e+O2(X,v=1)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(1)=data4me.rateExtra(i).value;
  case "e+O2(X,v=2)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(2)=data4me.rateExtra(i).value;
  case "e+O2(X,v=3)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(3)=data4me.rateExtra(i).value;
  case "e+O2(X,v=4)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(4)=data4me.rateExtra(i).value;
  case "e+O2(X,v=5)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(5)=data4me.rateExtra(i).value;
  case "e+O2(X,v=6)->e+O2(a1Dg),Excitation"
    k_O2Xi_O2a(6)=data4me.rateExtra(i).value;
  case "e+O2(X,v=1)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(1)=data4me.rateExtra(i).value;
  case "e+O2(X,v=2)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(2)=data4me.rateExtra(i).value;
  case "e+O2(X,v=3)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(3)=data4me.rateExtra(i).value;
  case "e+O2(X,v=4)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(4)=data4me.rateExtra(i).value;
  case "e+O2(X,v=5)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(5)=data4me.rateExtra(i).value;
  case "e+O2(X,v=6)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(6)=data4me.rateExtra(i).value;
  case "e+O2(X,v=7)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(7)=data4me.rateExtra(i).value;
  case "e+O2(X,v=8)->e+O2(b1Sg+),Excitation"
    k_O2Xi_O2b(8)=data4me.rateExtra(i).value;
  otherwise
    line=data4me.rateExtra(i).collDescription;
    if line(1:9)=='e+O2(X,v='
     num_lvl=line(10);
     num_sym=11;
     if line(11)~=')'
        num_lvl(2)=line(11);
        num_sym=12;
     end
     switch convertCharsToStrings(line(num_sym:end))
        case ")->e+2O(3P),Excitation"
            k_O2Xi_2O3P(str2num(num_lvl)+1)=...
                                data4me.rateExtra(i).value;
        case ")->O(-,gnd)+O(3P),Attachment"
            k_O2Xi__O3P_Om(str2num(num_lvl)+1)=...
                                data4me.rateExtra(i).value;
         otherwise
            if convertCharsToStrings(...
                                line(num_sym:num_sym+12))==")<->e+O2(X,v="
                num_lvl2=line(num_sym+13);
                if line(num_sym+14)~=')'
                    num_lvl2(2)=line(num_sym+14);
                end
                k_eV_f(str2num(num_lvl)+1, str2num(num_lvl2)+1)=...
                                            data4me.rateExtra(i).value(1);
                k_eV_b(str2num(num_lvl)+1, str2num(num_lvl2)+1)=...
                                            data4me.rateExtra(i).value(2);
                check_arr(str2num(num_lvl)+1, str2num(num_lvl2)+1)=...
                    check_arr(str2num(num_lvl)+1, str2num(num_lvl2)+1)+1;
            end
     end
    end
 end
end
% disp(num2str(k_eV))
% pause
    % e + O2(X, i) -> e + O2(X, i')
R(2:5)=(n1(1)*k_eV(:,1)-n1(2:5).*k_eV(:, 2))*ne;    % i=1-4
R(1)=-sum(R(2:5));                                  % i=0
    % e + O2(X,0) ->  e + 2O(3P)
R(6)=k_O2X_2O3P;
    % e + O2(X,0) <-> e + O2(a)
R(7)=ne*(n1(1)*k_O2X_O2a(1)-nO2a*k_O2X_O2a(2));
    % e + O2(X,0) <-> e + O2(b)
R(8)=ne*(n1(1)*k_O2X_O2b(1)-nO2b*k_O2X_O2b(2));
    % e + O2(a) <-> e + O2(b)
R(9)=ne*(nO2a*k_O2a_O2b(1)-nO2b*k_O2a_O2b(2));
    % e + O(3P) <-> e + O(1D)
R(10)=ne*(nO*k_O3P_O1D(1)-nO1D*k_O3P_O1D(2));
    % e + O2(X) -> e + O(3P) + O(1D)
R(11)=ne*(n1(1)*k_O2X__O3P_O1D(1));
    % e + O2(a1Dg) -> e + O(3P) + O(1D)
R(12)=ne*(nO2a*k_O2a__O3P_O1D(1));
    % e + O2(a) ->  e + 2O(3P)
R(13)=k_O2a__2O3P;
    % e + O2(b) ->  e + 2O(3P)
R(14)=k_O2b__2O3P;
    % e + O2(b) ->  e + O(3P) + O(1D)
R(15)=ne*nO2b*k_O2b__O3P_O1D;
    % e + O3(X) ->  e + O(3P) + O2(X,0)
R(16)=ne*nO3*k_O3X__O3P_O2X;
    % e + O2(X,0) -> 2e + O2(+,X)
R(17)=ne*n1(1)*k_O2X__e_O2p;
    % e + O2(a1Dg) -> 2e + O2(+,X)
R(18)=ne*nO2a*k_O2a__e_O2p;
    % e + O(3P) -> 2e + O(+,gnd)
R(19)=ne*nO*k_O3P__e_Op;
    % e + O2(X) -> 2e + O(3P) + O(+,gnd)
R(20)=ne*n1(1)*k_O2X__e_O3P_Op;
    % e + O2(a1Dg) -> 2e + O(3P) + O(+,gnd)
R(21)=ne*nO2a*k_O2a__e_O3P_Op;
    % e + O2(X) -> O(-,gnd) + O(3P)
R(22)=ne*n1(1)*k_e_O2X__Om_O3P;
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

% if mod(round(cputime), 10)==1
%     pr_vec(1:1+length(data4me.rate)+length(data4me.rateExtra))=0;
%     pr_vec(1)=t;
%     for i=1:length(data4me.rate)
%         pr_vec(1+i)=data4me.rate(i).value(1);
%     end
%     for i=1:length(data4me.rateExtra)
%         pr_vec(1+length(data4me.rate)+i)=data4me.rateExtra(i).value(1);
%     end
%     fprintf(fileKrate, [repmat('%d ', 1, numel(pr_vec)), '\n'], pr_vec);
% end

out=R';
% vibrational processes
    % e + O2(X, 1-6) -> e + O2(a)
R_O2Xi_O2a=ne*n1(2:7).*k_O2Xi_O2a';
    % e + O2(X, 1-8) -> e + O2(b)
R_O2Xi_O2b=ne*n1(2:9).*k_O2Xi_O2b';
    % e + O2(X, 1-41) -> e + 2O(3P)
k_O2Xi_2O3P(1)=0;
R_O2Xi_2O3P=ne*n1.*k_O2Xi_2O3P';
    % e + O2(X, 1-41) -> e + O(3P) + O(1D)
R_O2Xi__O3P_O1D=ne*n1.*k_O2X__O3P_O1D;
R_O2Xi__O3P_O1D(1)=0;
    % e + O2(X, 1-41) -> O(3P) + O-
R_O2Xi__O3P_Om=ne*n1.*k_O2Xi__O3P_Om';
R_O2Xi__O3P_Om(1)=0;
    % e + O2(X, 1-41) -> e + e + O2+
R_O2Xi_e_O2p=ne*n1*k_O2X__e_O2p;
R_O2Xi_e_O2p(1)=0;
    % eV: e + O2(X,v=i) <-> e + O2(X,v=j)
k_eV_f=k_eV_f.*n1;
k_eV_b=k_eV_f.*n1';
R_eV_cs=ne*(-sum(k_eV_f, 2)+sum(k_eV_b, 2)...
                                        +sum(k_eV_f, 1)'-sum(k_eV_b, 1)');
%     % MT test
% test=k_O2X_2O3P*Z;
% disp('test TM, first k rates')
% disp(k_O2X_2O3P)
% disp(test(1))
% disp('last k rates')
% disp(k_O2Xi_2O3P(end))
% disp(test(end))
% pause

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

function out=R_VV(M1, M2, T, n1, n2) % VV for O2-O2 collision
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

function out=R_VT(M, T, n1, na, config)

core=2:M.num_vibr_levels(1)-1; %c 1 по n-1 уровень

R=zeros(M.num_vibr_levels(1), 1);
% k_down_O2=k_VT(M, T, 1:M.num_vibr_levels-1, -1, M)';
kvt_temp=kvt_fho_di(T, 1, M);
k_down_O2=kvt_temp(1,:)';
k_down_O=kvt_temp(2,:)';
% k_up_O2=k_VT(M, T, 0:M.num_vibr_levels-2, 1, M)';
i=0:M.num_vibr_levels(1)-2;
delta_i=1;
dE=M.e_i(1,i+1)-M.e_i(1,i+1+delta_i);
k_up_O2=k_down_O2.*k_bf_VT(T, -dE)';
k_up_O=k_down_O.*k_bf_VT(T, -dE)';
% disp('krate FHO')
% disp(k_down_O2(1))
% disp(k_up_O2(1))

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

if config=="FHO10"
    for j=2:10
        kvt_temp=kvt_fho_di(T, j, M);
        k_down_O=kvt_temp(2,:)';
        R_temp=na*n1(1+j:end).*k_down_O;
        R(1+j:end)=R(1+j:end)-R_temp;
        R(1:end-j)=R(1:end-j)+R_temp;
        i=0:M.num_vibr_levels(1)-1-j;
        % delta_i=1;
        dE=M.e_i(1,i+1)-M.e_i(1,i+1+j);
        k_up_O=k_down_O.*k_bf_VT(T, -dE)';
        R_temp=na*n1(1:end-j).*k_up_O;
        R(1+j:end)=R(1+j:end)+R_temp;
        R(1:end-j)=R(1:end-j)-R_temp;
    end
end

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