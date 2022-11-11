k = 1.3807e-23; % постоянная Больцмана (Дж/К)
N_a=6.02214076e23;  % Avogadro constant

    % Oxygen atoms
O.mass=2.65676E-26;
O.form_e = 4.098045681049634e-19;
O.num_elex_levels=1;
O.e_E=[0 3.15203E-19];    % O(3P), O(1D)

    % oxygen molecules
O2_O2_coll_diam=3.458e-10; % диаметр столкновения O2 + O2
O2.num_elex_levels=3;       % number of elecronical levels 
                            %   (O2(X), O2(a), O2(b))
O2.num_vibr_levels=[42 1 1];  % number of vibrational levels
% O2.num_vibr_levels=[1 1 1];   % test
O2.mass=5.31353E-26;
O2.diameter=3.458e-10;
O2.sigma=2;                 % symmetry factor
O2.Be=[143.768 142.64 140.03699999999998];
O2.mltpl_atoms_mass=O.mass*O.mass;
O2.s_e=[3 2 1];               % electronic statistical weight
% O2.we=0;                    % spectroscopic const 
    % electron excited levels energy from Savelev tables
O2.e_E=[0 1.57289E-19 2.62114E-19]; 
O2.red_osc_mass=1.13923e-26;
O2.mA=0.5;  O2.mB=0.5;
O2.form_e=0;
% O2.diss_e=59379*k;
O2.diss_e=[8.19609E-19 8.20325E-19 8.35261E-19];
E_vibr_levels_eV=[  % колебательная энергия в eV. 41 уровень.
0 0.000 14 2.435 28 4.280
1 0.196 15 2.587 29 4.382
2 0.388 16 2.735 30 4.476
3 0.573 17 2.881 31 4.565
4 0.756 18 3.024 32 4.651
5 0.937 19 3.164 33 4.730
6 1.117 20 3.301 34 4.794
7 1.291 21 3.436 35 4.847
8 1.461 22 3.568 36 4.898
9 1.629 23 3.696 37 4.938
10 1.796 24 3.821 38 4.960
11 1.960 25 3.942 39 4.976
12 2.122 26 4.059 40 4.987
13 2.281 27 4.172 41 4.994];
e_vibr=zeros(1, O2.num_vibr_levels(1));
e_vibr(1:14)=E_vibr_levels_eV(:,2);
e_vibr(15:28)=E_vibr_levels_eV(:,4);
e_vibr(29:42)=E_vibr_levels_eV(:,6);
e_vibr=e_vibr*1.602176487e-19;   % колебательная энергия в Дж
O2.e_i=e_vibr;
O2.e_0=1.559405589690000e-20;   %предположительное значение энергии 
                                %нулевого уровня, взятое из ангарм. осц.
O2.e_i(2, :)=0;
O2.e_i(3, :)=0;
    % ozon moleculas
O3.mass=3*15.999*1.660539040e-27;
O3.e_i=1553*1.38064878066922E-23;

    % electron data
e.mass=9.1093837015e-31;

    % overall
Prcl.O2=O2; Prcl.O=O; Prcl.O3=O3; Prcl.e=e;

%%%%%%%%%%%%%%%% COLLISIONS %%%%%%%%%%%%%%
    % O2-O2
Coll_O2_O2.ArrA=2e21/N_a*1e-6;    Coll_O2_O2.ArrN=-1.5;
    % O2-O
Coll_O2_O.ArrA=1e22/N_a*1e-6;     Coll_O2_O.ArrN=-1.5;
    % overall
Coll.O2_O2=Coll_O2_O2;  Coll.O2_O=Coll_O2_O;