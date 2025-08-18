function k_down = kvt_Esposito(T)%, M1, M2, ind_e, anharm)
% Esposito k_VT
% only dv = 1
% O2-O
% 15.08.2025 by Maksim Melnik

% load("..\data\models_data\Esposito_VT_O2O.mat", "b_d1_arr")

%% b_ij
temp_b1 = [
-26.18993227    
-1.69828917     
3.349076e19     
-3.946126e20    
1.391056e19     
];

% j,k Dv = 1 Dv = 2-10 Dv = 11-20 Dv = 21-30
temp_b2 = [
7.83331061      
3.71221451      
3.573261e20     
6.433503e20     
-2.901352e19    
];

% j,k Dv = 1 Dv = 2–10 Dv = 11–20 Dv = 21–30
temp_b3 = [
0.37163948      
0.10587091      
-5.312491e19    
3.754092e19     
-1.189832e18    
];

b_d1_arr = [temp_b1'; temp_b2'; temp_b3'];

% dv = 1;
v = 1:46;
DegF = 1/27;

% a = zeros(3, length(v));
a = b_d1_arr(:, 1) + b_d1_arr(:, 2) .* log(v) + ...
    (b_d1_arr(:, 3) + b_d1_arr(:, 4) .* v + b_d1_arr(:, 5) .* v .* v) ...
                                                    ./ (1e21 + exp(v));
k_down = DegF * exp(a(1, :) + a(2, v)/log(T) + a(3, v)*log(T));
end