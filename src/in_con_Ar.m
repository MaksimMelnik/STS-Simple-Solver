function Y = in_con_Ar(IN)
% format long e
V_K = 1.380649e-23;
syms n v t % dim-less variables

% conditions in free stream
mCO = IN(1); % mass of CO molecule, kg
v0 = IN(2); % velocity, m/s
T0 = IN(3); % temperature, K
mAr = IN(4); % mass of Ar atom, kg
f = IN(5); % fraction of CO in the mixture

% flow parameters just behind a shock
F=mAr*(1-f)+mCO*f;
C2 = F*v0^2/V_K/T0;
% 
S = [n*v == 1,...
     n*t+n*v^2*C2 == 1+C2,...
     (2.5+f)*t+v^2*0.5*C2 == 2.5+f+0.5*C2];

N = vpasolve(S,[n,t,v],[10,100,0.1]);
X = [N.n;N.t;N.v];
Y = double(X);
gg=[2 4 6];
Z=(Y-1).^2;
if sum(Z(gg-1))>sum(Z(gg))
    gg=gg-1;
end
Y=Y(gg);
% Y = Y([2 4 6]);
% Y = Y([1 3 5]);