clearvars; close all; clc;
addpath('../data');

% Количество запусков модели (вообще 10, но можно поэкспериментировать)
num_runs_sens = 10;


%% выбор лучшего LHS

num_candidates = 5000;   % сколько разных LHS генерировать
d = 3;                   % размерность параметров
n = num_runs_sens;      

fprintf('=== Поиск лучшего LHS для устойчивой метамодели ===\n');

lP_0 = @(x) 1;
lP_1 = @(x) sqrt(3) * x;
lP_2 = @(x) sqrt(5)/2 * (3*x.^2 - 1);
best_score   = -Inf;
best_LHS  = [];
best_det  = -Inf;
best_cond = Inf;

for k = 1:num_candidates
    
    % ---------- 1. Генерируем LHS ----------
    X = lhsdesign(n, d, 'criterion','maximin');  % в [0,1]
    X = 2 * X - 1;   % масштабируем в [-1,1]
  
    Psi = zeros(n, 10);
    
    for i = 1:n
        xi = X(i,:);
        
    Psi(i, 1) = lP_0(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));      
    Psi(i, 2) = lP_1(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi(i, 3) = lP_0(xi(1)) * lP_1(xi(2)) * lP_0(xi(3)); 
    Psi(i, 4) = lP_0(xi(1)) * lP_0(xi(2)) * lP_1(xi(3)); 
    Psi(i, 5) = lP_2(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi(i, 6) = lP_0(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));  
    Psi(i, 7) = lP_0(xi(1)) * lP_0(xi(2)) * lP_2(xi(3)); 
    Psi(i, 8) = lP_1(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));  
    Psi(i, 9) = lP_1(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));  
    Psi(i, 10) = lP_0(xi(1)) * lP_1(xi(2)) * lP_1(xi(3)); 
        end
    
    G = Psi' * Psi; % информационная матрица

    detG  = det(G);
    condG = cond(G);
    
    if condG < 100         
        score = detG;      % максимизируем определитель
        if score > best_score
            best_score = score;
            best_LHS   = X;
            best_det   = detG;
            best_cond  = condG;
        end
    end
end

fprintf('Лучший LHS найден!\n');
fprintf('det(Psi^T Psi) = %e\n', best_det);
fprintf('cond(Psi^T Psi) = %e\n', best_cond);

LHS_points = best_LHS;

%% основной код

results_matrix = zeros(501, 11);  

for run_idx = 1:num_runs_sens
    fprintf('\n=== ЗАПУСК %d из %d ===\n', run_idx, num_runs_sens);
    
    current_point = LHS_points(run_idx, :);
    factors = 10.^(2 * current_point);
    
    % Загружаем исходные данные
    load('../data/particles.mat', 'O2', 'O', 'Ar');
    
    % Меняем параметры
    O2.diss_Arrhenius_A('O')  = O2.diss_Arrhenius_A('O')  * factors(1);
    O2.diss_Arrhenius_A('O2') = O2.diss_Arrhenius_A('O2') * factors(2);
    O2.diss_Arrhenius_A('Ar') = O2.diss_Arrhenius_A('Ar') * factors(3);
    save 'particles1.mat' O2 O Ar;
    
    fprintf('Параметры запуска %d:\n', run_idx);
    fprintf('  O2+O:  %e\n', O2.diss_Arrhenius_A('O'));
    fprintf('  O2+O2: %e\n', O2.diss_Arrhenius_A('O2'));
    fprintf('  O2+Ar: %e\n', O2.diss_Arrhenius_A('Ar'));
    
    % Сохраняем состояние
    save('temp_sensitivity_vars.mat', 'run_idx', 'results_matrix', 'num_runs_sens', 'LHS_points');
    
    % Запускаем черный ящик
    fprintf('Запускаем O2_Ar_SW_Streicher...\n');
    O2_Ar_SW_Streicher;
    
    % Загружаем обратно
    load('temp_sensitivity_vars.mat');
    
    
    time = time_ms_1;
    concentration_O2 = n_O2_1 / Na;
    
    if run_idx == 1
        results_matrix(:, 1) = time;
    end
    results_matrix(:, run_idx + 1) = concentration_O2;
    
    fprintf('Запуск %d завершен\n', run_idx);
end


if exist('temp_sensitivity_vars.mat', 'file')
    delete('temp_sensitivity_vars.mat');
end

% Сохраняем результаты
save('sensitivity_results_LHS.mat', 'results_matrix');
writematrix(results_matrix, 'sensitivity_results_LHS.txt', 'Delimiter','tab');

%%
%%
figure;
plot(results_matrix(:,1), results_matrix(:,2:end));

%% решение слу
lP_0 = @(x)  legendreP(0,x);
lP_1 = @(x) sqrt(3) *  legendreP(1,x);
lP_2 = @(x) sqrt(5) *  legendreP(2,x);
Np = 9;
Nt = Np + 1;

Psi = zeros(num_runs_sens, Nt);

for i = 1:num_runs_sens
    xi = LHS_points(i,:);
    Psi(i, 1) = lP_0(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));      
    Psi(i, 2) = lP_1(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi(i, 3) = lP_0(xi(1)) * lP_1(xi(2)) * lP_0(xi(3)); 
    Psi(i, 4) = lP_0(xi(1)) * lP_0(xi(2)) * lP_1(xi(3)); 
    Psi(i, 5) = lP_2(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi(i, 6) = lP_0(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));  
    Psi(i, 7) = lP_0(xi(1)) * lP_0(xi(2)) * lP_2(xi(3)); 
    Psi(i, 8) = lP_1(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));  
    Psi(i, 9) = lP_1(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));  
    Psi(i, 10) = lP_0(xi(1)) * lP_1(xi(2)) * lP_1(xi(3)); 
end

alpha = zeros(Nt, 501);

for i = 1:501
    R = results_matrix(i, 2:end)';
    alpha(:,i) = linsolve(Psi, R);
end

%% индексы Соболя

IN_Sobol = zeros(500, 4);
IN_Sobol(:,1) = time(2:501);

for i=2:501
    D = sum(alpha(2:end, i).^2);
    IN_Sobol(i-1,2) = (alpha(2,i)^2 + alpha(5,i)^2 + alpha(8,i)^2 + alpha(9,i)^2)   / D;
    IN_Sobol(i-1,3) = (alpha(3,i)^2 + alpha(6,i)^2 + alpha(8,i)^2 + alpha(10,i)^2) / D;
    IN_Sobol(i-1,4) = (alpha(4,i)^2 + alpha(7,i)^2 + alpha(9,i)^2 + alpha(10,i)^2)/ D;
end

figure;
plot(IN_Sobol(:,1), IN_Sobol(:,2:4));
legend('A_O', 'A_O2', 'A_Ar');
xlabel('t');
ylabel('Total Sobol Indices');


