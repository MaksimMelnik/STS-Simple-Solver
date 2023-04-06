function Y = in_con_O2(IN)
% Solving initial conditions behind shock wave.
format long e
k=1.380649e-23;
syms n v t % dim-less variables

% conditions in free stream
m = IN(1); % mass of molecule, kg
v0 = IN(2); % velocity, m/s
T0 = IN(3); % temperature, K

% flow parameters just behind a shock
C1 = m;
C2 = C1*v0^2/k/T0;
% 
S = [n*v == 1,...
     n*t+n*v^2*C2 == 1+C2,...
     3.5*t+v^2*0.5*C2 == 3.5+0.5*C2];

N = vpasolve(S,[n,t,v],[10,100,0.1]);
X = [N.n;N.t;N.v];
Y = double(X);
gg=[2 4 6];
Z=(Y-1).^2;
if sum(Z(gg-1))>sum(Z(gg))
    gg=gg-1;
end
Y=Y(gg);
%Y = Y([2 4 6]);
%Y = Y([1 3 5]);