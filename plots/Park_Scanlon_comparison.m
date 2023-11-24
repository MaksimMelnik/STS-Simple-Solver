Na = 6.02214076e23; 
k = 1.380649e-23;
reaction=["O2 + N → O + O + N"; "O2 + NO → O + O + NO"; ...
    "O2 + N2 → O + O + N2"; "O2 + O2 → O + O + O2";...
    "O2 + O → O + O + O"; "N2 + O → N + N + O";...
    "N2 + O2 → N + N + O2"; "N2 + NO → N + N + NO";...
    "N2 + N2 → N + N + N2"; "N2 + N → N + N + N";...
    "NO + N2 → N + O + N2"; "NO + O2 → N + O + O2";
    "NO + NO → N + O + NO"; "NO + O → N + O + O";
    "NO + N → N + O + N"];
preExp_Park = [1e22; 2e21; 2e21; 2e21; 1e22; 3e22; 7e21; 7e21; 7e21; 3e22;...
    5e15; 5e15; 1.1e17; 1.1e17; 1.1e17]/Na*1e-6;
preExp_Scanlon = [1.1e-10; 1.1e-10; 1.3e-10; 5.33e-11; 1.5e-10; 4e-12; 1.5e-11;...
    1.5e-11; 4.1e-12; 1e-11; 2.1e-10; 2e-10; 1e-10; 4e-10; 4e-10];
n_Park=[-1.5; -1.5; -1.5; -1.5; -1.5; -1.6; -1.6; -1.6; -1.6; -1.6;...
    0; 0; 0; 0; 0];
n_Scanlon=[-1; -1; -1; -1; -1.05; -0.54; -0.68; -0.68; -0.62; -0.68; -1;...
    -1; -1; -1.1; -1.1];
Ea_Park =  [59750; 59750; 59750; 59750; 59750; 113200; 113200; 113200;...
    113200; 113200; 75500; 75500; 75500; 75500; 75500]*k;
Ea_Scanlon = [8.197; 8.197; 8.197; 8.197; 8.197; 15.67; 15.67; 15.67; ...
    15.67; 15.67; 10.43; 10.43; 10.43; 10.43; 10.43]*1e-19;

T=2000:10:14000;
kf_Park=preExp_Park.*T.^n_Park.*exp(-Ea_Park./T/k);
kf_Scanlon=preExp_Scanlon.*T.^n_Scanlon.*exp(-Ea_Scanlon./T/k);
const=table(reaction, preExp_Park, preExp_Scanlon, n_Park, n_Scanlon, Ea_Park, Ea_Scanlon)

for i=1:15
figure
semilogy(T, kf_Park(i,:), 'k-');
hold on
semilogy(T, kf_Scanlon(i,:), 'k--');
hold off
title("Reaction: " + reaction(i));
xlabel("T, K");
ylabel("log K_f");
legend("Park model", "Scanlon model", 'Location','best');
exportgraphics(gca,reaction(i)+".jpg");
end
