% ��� ��� ������� 0d ���������� ������������������ ������ � �������.
% �������� ��������� ����-���� �� ����� � ������ ������������� ��������.
% �� ������ code_Tiago_imitation, �� ������ � �������������� ��������
% 12.02.2021

% ������:

% ��������� FHO ��� ����������
% �������� ����������� �� ������� ������ (��������� ������ �����)
    % 1. O2(X,41)+O2(X)->2O(3P)+O2(X)
    % 2. O2-O2 VV �����������
    % 3. ����������� �� ������ O(3P) � O(1D)
% ������ �� �����, ������ VT O ���� ������, ��� � ������
% �������� ���������� ����������� ��������
    % �������� � ������������ ��� �� �����
    % �� ����� ������� �� ���� �������� ��������� ������������ �����
    %   ������������ � �� ��� ��� �������� ��������?

% �������� ���� �� ������-������ ������ (��� ������������ � O2a � ������
%   ������)
% ���������� ��� �������
    % ������ ����� ����
% ������������� � ������������� �������� ����
% ���������� FHO
% ����� ����������� ������������ ��� ���������, ���������� � 
%   ������������� ������� (������ ����� (��� ������ ��� ��������), ����� 
%   ����� � ������ ��� � ��������� �������)
    % 1. e + O2(X, 1-41) -> e + O(3P) + O(3P)
    % 2. e + O2(X, 1�41) -> e + O(3P) + O(1D)
    % 3. e + O2(X, 1-41) -> e + e + O2 +
    % 4. e + O2(X, 1-6) -> e + O2(a)
    % 5. e + O2(X, 1-8) -> e + O2(b)
    % 6. e + O2(X, 1-41) -> O(3P) + O-
% ����������� ��������� ��������������� ��� ������������� ������ ������
% ��������, ��� ���������� �� ������ (8, 8.3 � ������ ������)
% ����������� � �������������: �������� � ����� ��� ��� �� �
% O(-,gnd) + O(3P) -> e + O2(X)     ��� krate �� ��������� � ������ �
%   � ����� (����� ������� ���, ���������)
% � ��������� ��������� � ������� ������� -- �2(�����), � � ����� ������
%   �2(�). ���������� � ��, ����� �� ������� ���.
% ��������� �������� ������ ��� � ���� ����� ����
% ����������� ������� ���� ����������� ��������� �������
% ����������� �������� ������ ��������������
% ��������, ��� ��������� ����������� ����� ����������� � ���������
%   �����������
% ��������� ����� � ����, ������� ����� ��������� ������� ��� ������
%   ������, ����� �� ����������� Qin
% ���������� Qin ��� ������ ���������
% ����������� e+O2(a)->e+2O(3P) �� ��������, � �� ��� ������
% �������� Qin ��� ���������, �� ��������� � �����������
%   (� ����� ����������� �� ��������� ���)
    % �������� ��������
    % ����������
    % �������� ������
    % �������� � �� ��� ����������� (����)
% ���������� ��������� ���������� �� ��������� ������
% �������� ������������� ��������
% ������ ����������� � ������������ �� ������, ������ ��� ����� �������
%   ������ ���� � ������������� �� ����
    % ������ � -- ���� ������ ������ ������ � ���������, �� � ������� 
    %   ��������, ����� ��� ������ � ������
    % ������ ����� ������� ������, ���������, ����� �� ��������� 
    %   ������������� ��������� � ��������� ������� ��� ������������
% ��� �������� � ������� ����
% ��� � ��� ����������� ����� ���
    

fclose('all');
checkX=1e-20; %����� �� ������� ����� �������� ������� � �������
ct=cputime-10;

%   ���������
k = 1.3807e-23; % ���������� ��������� (��/�)
R = 8.3144621;    % ������������� ������� ����������
h = 6.626070041e-34;    % ���������� ������ (��/�)
HBAR = 1.054571800e-34;
c = 299792458;
Torr = 133.332;%.322368; % ����� �� � 1 ���
alpha_FHO =4e10;
E_FHO =200;
SVT_FHO=0.44;

    % ��������� ��� ��������� ���� ������ ���� � ������������ ����� ����
species_ini_STS

            
          % p,Torr  fO2  fO2* fO  T   EdN Ion. d.
incon_mas=[ 1       1    0    0   400 50  0
            1       0.7  0.1  0.2 400 50  0
            10      1    0    0   400 50  0
            10      0.7  0.1  0.2 400 50  0
            0.01    1    0    0   620 148 9.7e-5 % first experiment case
            0.01    1    0    0   700 116 1e-4   % second
            0.01    1    0    0   780 96  6.6e-5 % third
            0.01    1    0    0   900 83  2.9e-5 % forth
            0.01    1    0    0   400 50  0      % Polyh. tests
        ];
test_ar_p=[0.01 0.1 1 10];
test_ar_T=[300 650 1000];
ip=3; iT=3;
% incon_mas(9, :)=[test_ar_p(ip) 0.989 0.001 0.01 test_ar_T(iT) 70 0];
incon_mas(9, :)=[test_ar_p(ip) 1 0 0 test_ar_T(iT) 70 0];
ind=3;  % � �������� ���������� ����� �� ����������
T=incon_mas(ind, 5);
p=incon_mas(ind, 1)*Torr;     % �������� � ��������
EdN=incon_mas(ind, 6);        % ����
n_all=p/(k*T);                % ���������� ���� �������
n1=n_all*incon_mas(ind, 2);   % ���������� ������� O2(X)
nO2a=n_all*incon_mas(ind, 3); % ���������� ������� O2(a)
nO=n_all*incon_mas(ind, 4);   % �������� ��������� ������ ��������� O

ne=6.5e15;                  % �������� ��������� ���������� � �-3
is_e=1;
if is_e==0
    ne=0;
end
Tv=T;                       % ������������� �����������
T0=T;                       % ����������� ��� ����������������
n0=n1;                      % �� �� ����� � n0 ��� ����������������
tau=k/(p*pi*O2_O2_coll_diam*O2_O2_coll_diam)*sqrt(32*T/(6*R));

% ������������� ���������
	Zv=sum(exp(-O2.e_i(1,:)/Tv/k));
    n=n1/Zv.*exp(-O2.e_i(1,:)/Tv/k);
if ind<5
    n=n*0;
    n(1)=n1;
end
    
E0=1.5*(n1+nO)*k*T+n1*k*T+(O2.e_i(1,:)+O2.e_0)*n'+O.form_e*nO;

num_O2=1;%O2.num_vibr_levels(1);
y=n/n0;                         % ������� ������������� ������
num_O2a=num_O2;
if O2.num_elex_levels>1
    num_O2a=num_O2+O2.num_vibr_levels(1);
    y(num_O2a)=n_all*incon_mas(ind, 3)/n0; % first electron excited state
end
num_O2b=num_O2a;
if O2.num_elex_levels>2
    num_O2b=num_O2a+O2.num_vibr_levels(2);
    y(num_O2b)=0/n0;                       % second electron excited state
end
num_O2p=num_O2b+O2.num_vibr_levels(3);
num_O=num_O2p+1;
num_O1D=num_O+1;
num_Op=num_O1D+1;
num_Om=num_Op+1;
num_O3=num_Om+1;
num_O3exc=num_O3+1;
num_e=num_O3exc+1;
num_T=num_e+1;
num_E_e=num_T+1;

y(num_O2p)=0;%n1*incon_mas(ind, 7)/n0; % ion O2+
y(num_O)=nO/n0;     % �����
y(num_O1D)=0;
y(num_Op)=0;        % ion O+
y(num_Om)=0;        % ion O-
y(num_O3)=0;
y(num_O3exc)=0;     % O3 excited
y(num_e)=ne/n0;     % ���������
y(num_T)=T/T0;
% y(num_E_e)=0;%Tv/T0;

% �������� ��������
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-70);
% options = odeset('RelTol', 3e-14, 'AbsTol', 1e-60);     % ��� ������������
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-50);
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-40);
% options = odeset('RelTol', 3e-14, 'AbsTol', 1e-30);
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-20);     % �����
% options = odeset('RelTol', 3e-14, 'AbsTol', 1e-15);
% options = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
% options = odeset('RelTol', 3e-13, 'AbsTol', 1e-13);
options = odeset('RelTol', 3e-8, 'AbsTol', 1e-8);
if is_e==0
    options = odeset('RelTol', 3e-14, 'AbsTol', 1e-100);
    options = odeset('RelTol', 3e-14, 'AbsTol', 1e-80);
%     options = odeset('RelTol', 3e-14, 'AbsTol', 1e-40);
%     options = odeset('RelTol', 3e-14, 'AbsTol', 1e-15);
%     options = odeset('RelTol', 1e-8, 'AbsTol', 1e-15);
%     options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
end
% xspan_m=[0  1e2
%          0  1e2
%          0  1e2    %5.4e-8 1e-6 % 3e-1
%          0  1e2
%          0  5e0];
 xspan_m=[0  1e1];  % 1e2
%  if iT==3
%      xspan_m=[0  1e5];  % 1e2
%  end
% xspan=xspan_m(ind, :)/tau;
xspan=xspan_m/tau;
% disp(xspan(2)*tau)
% pause
fileData = fopen('print_data_test.txt','w');
fileKrate = fopen('print_krate.txt','w');
fprintf(fileData, ['%% t     nO2(X, 0-41)    nO2(a)  nO2(b)  nO2+'...
                    '   nO(3P)  nO(1D) nO+  nO- O3  O3exc   ne  T\n']);
% load LoKI-master\Code\Output\RNF\data4Maksim.mat
% data4me=data4Maksim;
% pr_vec='%t';
% for i=1:length(data4me.rate)
%     pr_vec=[pr_vec '\t' data4me.rate(i).collDescription];
% end
% for i=1:length(data4me.rateExtra)
%     pr_vec=[pr_vec '\t' data4me.rateExtra(i).collDescription];
% end
% fprintf(fileKrate, [pr_vec '\n']);
if is_e==1
run_LoKI_STS(T, n'/n0, y(num_O2a), y(num_O2b), ne/n0, nO/n0, ...
                                                y(num_O1D), n0, EdN, 0);
end

config.V="Paper";
config.V="FHO";
% config.V="FHO10";
config.LoKI=0;  % -1 no LoKI, 0 every step, N>0 -- how often to start LoKI
%{
% param.VTO2O="Esposito";
% param.VTO2O="FHO";
% param.VTO2O="FHO10";
% param.VTO2O2="Annusova";
% param.VTO2O2="FHO";
% param.VTO2O2="FHO10";
% param.VVO2O2="Annusova";
% param.VVO2O2="FHO";
% param.VVO2O2="FHO10";
%}

disp('ODE solver started.')
disp(['T=' num2str(T) ', Tv=' num2str(Tv) ', EdN=' num2str(EdN)...
    ', p=' num2str(p/Torr) ...
    ', Electrons on: ' num2str(is_e) ','])
disp(['ind=' num2str(ind) ', n1/n0=' num2str(n(2)/n(1)) ...
    ', xspan=[' num2str(xspan(1)*tau) ', ' num2str(xspan(2)*tau) ']'])
disp(config.V)
% disp([y(num_T)])% y(num_E_e)])
tic
[X, Y]= ode15s(@(t,y) Rpart_Tiago_STS(t, y, Prcl, Coll, tau, ...
    xspan, T0, n0, is_e, EdN, fileData, fileKrate, config),...
    xspan, y, options);             % ������������� ������������ ode15s
X=X*tau;
Y(:, 1:num_e)=Y(:, 1:num_e)*n0;
Y(:, num_T)=Y(:, num_T)*T0;
data_Tiago_test=[X Y];

% run_LoKI(T, Y(end,num_O2)/n0, Y(end, num_O2a)/n0, Y(end, num_O2b)/n0, ...
%     ne/n0, Y(end, num_O)/n0, ...
%                                     Y(end, num_O1D)/n0, n0, EdN, 1e2);

toc
fclose('all');
disp(['p=' num2str(p/Torr) ', T=' num2str(T) ', Tv0=' num2str(Tv) ...
    ', n1/n0=' num2str(n(2)/n(1))])
disp(['ne=' num2str(Y(end, num_e), '%.4e')])
disp(config.V)
%%
% temp=data_Tiago_test;
e_i=O2.e_i(1,:);
num_vibr_levels=O2.num_vibr_levels(1);
% zero_p=sum(temp(1, 2:num_vibr_levels+1))*temp(1,end);
n_o2=sum(Y(:, 1:num_O2p), 2);
n_o=sum(Y(:, num_O:num_Om), 2);
n_o3=sum(Y(:, num_O3:num_O3exc), 2);
% all_p=n_o2.*temp(:,end);
rho=n_o2*O2.mass+n_o*O.mass+n_o3*O3.mass;
disp(['max error rho ' num2str(max(abs((rho(1)-rho(:))/rho(1))))])
% tr vibr form rot
T_vec=Y(:, num_T);
% Tv = O2.e_i(1,2)./(k*log(temp(:, 2)./temp(:, 3)));
Energy=1.5*(n_o2+n_o+n_o3)*k.*T_vec+...
    sum(([O2.e_i(1,1:O2.num_vibr_levels(1)) ...
    O2.e_i(2, 1:O2.num_vibr_levels(2)) ...
    O2.e_i(3,1:O2.num_vibr_levels(3))...
    ]+O2.e_0).*Y(:,1:num_O2b+O2.num_vibr_levels(3)-1),2)...
      + n_o*O.form_e+(n_o2+n_o3)*k.*T_vec+Y(:, num_O3exc)*O3.e_i;
disp(['max error Energy ' ...
    num2str(max(abs((Energy(1)-Energy(:))/Energy(1))))])
%% added processes:
    % ��� nm ����, �� �� �������� �� ������ nO2
    % ������������� �������� ��� �� ���������
%   O2(X)
% O2(X,i)   +O2(X,j)    ->O2(X,i-1) +O2(X,j+1) % and reverse | vibr, VV
% O2(X,i)   +O2(X)      ->O2(X,i')  +O2(X)     % vibr, VT
% O2(X,i)   +O(3P)      ->O2(X,i-di)+O(3P)     %di=1:30+reverse | vibr, VT
% O2(X,i)   +O2(X)      ->2O(3P)    +O2(X)     % diss
% % O2(X)    +O2(a)      ->2O(3P)    +O2(a)    % diss
% % O2(X)    +O2(b)      ->2O(3P)    +O2(b)    % diss
% O2(X,i)   +O(3P)      ->2O(3P)    +O(3P)     % diss
% O2(a)     +O2(X)      ->2O(3P)    +O2(X)     % diss
% O2(b)     +O2(X)      ->2O(3P)    +O2(X)     % diss
% 2O(3P)    +O2(X)      ->O2(X)     +O2(X)     % rec
% % 2O(3P)   +O2(a)      ->O2(X)     +O2(a)    % rec
% % 2O(3P)   +O2(b)      ->O2(X)     +O2(b)    % rec
% 2O(3P)    +O2(X)      ->O2(a)     +O2(X)    % rec
% 2O(3P)    +O2(X)      ->O2(b)     +O2(X)    % rec
% 2O(3P)    +O(3P)      ->O2(X)     +O(3P)    % rec
% O2(b)     +O(3P)      ->O2(X)     +O(3P)    % deex
% O2(a)     +O2(X)      ->O2(X)     +O2(X)    % deex
% % O2(a)    +O2(a)      ->O2(X)    +O2(a)    % deex
% % O2(a)    +O2(b)      ->O2(X)    +O2(b)    % deex
% O(1D)     +O2(X)      ->O(3P)     +O2(a)    % mix (deex, exci)
% O(1D)     +O2(X)      ->O(3P)     +O2(b)    % mix (deex, exci)
% O(1D)     +O2(X)      ->O(3P)     +O2(X)    % deex
% O2(a)     +O2(a)      ->O2(b)     +O2(X)    % exch
% O2(a)     +O(3P)      ->O2(X)     +O(3P)    % deex
% O(1D)     +O3(X)      ->2O2(X)              % exch (low)
% O(1D)     +O3(X)      ->O2(X)     +2O(3P)   % exch (low)
% O(3P)     +O2(X)+O(3P)->O3(X)     +O(3P)    % rec
% O(3P)     +O2(X)+O2(X)->O3(X)     +O2(X)    % rec
% O(3P)     +O2(X)+O3(X)->O3(X)     +O3(X)    % rec (low)
% O2(a)     +O3(X)      ->2O2(X)    +O(3P)    % mix (diss, deex)
% O2(b)     +O3(X)      ->2O2(X)    +O(3P)    % mix (diss, deex)
% O(3P)     +O3(X)      ->2O2(X)              % exch
% O(3P)     +O3(X)      ->O2(a)     +O2(X)    % exch
% O(3P)     +O3(X)      ->O2(b)     +O2(X)    % exch
% O(+)      +O3(X)      ->O2(+,X)   +O2(X)    % exch
% O(+)      +O2(X)      ->O2(+,X)   +O(3P)    % exch
% O(-)      +O(3P)      ->O2(X)     + e       % mix (rec, deion)
% O(-)      +O2(X)      ->O3(X)     + e       % mix (rec, deion)
% O(-)      +O2(b)      ->O(3P)     +O2(X)+e  % mix (deion, deexci)
% O2(+,X)   +O(-)       ->O2(X)     +O(3P)    % deion
% O(3P)     +O2(X)+O2(X)->O3(exc)   +O2(X)    % rec
% O3(exc)   +O2(X)      ->O3(X)     +O2(X)    % deex
% O2(a)     +O3(exc)    ->2O2(X)    +O(3P)    % mix (diss, deex)
% O(3P)     +O3(exc)    ->2O2(X)              % exch
% O2(X)     + e         ->2O(3P)    + e       % elec, diss
% O2(X)     + e        <->O2(a)     + e       % elec, ex
% O2(X)     + e        <->O2(b)     + e       % elec, ex
% O2(X)     + e         ->O(3P)     +O(1D)+e  % elec, diss
% O3(X)     + e         ->O(3P)     +O2(X)+e  % elec, diss
% O2(X)     + e         ->O2(+,X)   + 2e      % elec, ion
% O2(X)     + e         ->O(3P)     + O(+)+2e % elec, mix (diss, ion)
% O2(X,1-6) + e         ->O2(a)     + e       % elec, vibr, ex
% O2(X,1-8) + e         ->O2(b)     + e       % elec, vibr, ex
% O2(X,1-41)+ e         ->2O(3P)    + e       % elec, vibr, diss
% O2(X,1-41)+ e         ->O(3P)     +O2(1D)+e % elec, vibr, diss
% O2(X,1-41)+ e         ->O2+       + 2e      % elec, vibr, ion
% O2(X,1-41)+e          ->O(3P)     + O-      % elec, vibr, ion
% 	O2(a)
% % O2(a) +O2(a)->2O(3P) +O2(a)    % diss
% % O2(a) +O2(b)->2O(3P) +O2(b)    % diss
% O2(a)     +O(3P)->2O(3P) +O(3P)    % diss 
% % 2O(3P)+O2(a)->O2(a)  +O2(a)    % rec 
% % 2O(3P)+O2(b)->O2(a)  +O2(b)    % rec
% % 2O(3P)+O2(a)->O2(b)  +O2(a)    % rec
% O2(b)     +O(3P)->O2(a)  +O(3P)    % deex
% O(+)      +O2(a)->O2(+,X)+O(3P)    % exch
% O(-)      +O2(a)->O3(X)  + e       % mix (rec, ion)
% O2(a)     + e   ->2O(3P) + e       % elec, diss
% O2(a)     + e  <->O2(b)  + e       % elec, ex
% O2(a)     + e   ->O(3P)  +O(1D)+e  % elec, diss
% O2(a)     + e   ->O2(+,X)+ 2e      % elec, ion
% O2(a)     + e   ->O(3P)  +O(+) +2e % elec, mix (diss, ion)
% O2(X)     + e   ->O(-)   +O(3P)    % elec, mix (diss, ion)
% O2(a)     + e   ->O(-)   +O(3P)    % elec, mix (diss, ion)
%   O2(b)
% % O2(b) +O2(b)->2O(3P)+O2(b)    % diss
% O2(b) +O(3P)->2O(3P)+O(3P)    % diss
% % 2O(3P)+O2(b)->O2(b) +O2(b)    % rec
% O2(b) + e   ->2O(3P)+ e       % elec, diss
% O2(b) + e   ->O(3P) +O(1D)+e  % elec, diss
%   O(3P)
% O(1D)  +O(3P)->O(3P) + O(3P)   % deex
% O(+)   +O(-) ->2O(3P)          % deion
% O3(exc)+O(3P)->O3(X) + O(3P)   % deex
% O(3P)  + e  <->O(1D) +e        % elex, exci
% O2(+,X)+ e   ->2O(3P)          % elec, diss
% O2(+,X)+ e   ->O(3P) + O(1D)   % elec, diss
% O(3P)  + e   ->O(+)  + 2e      % elec, ion
% O(-)   + e   ->O(3P) + 2e      % elec, deion