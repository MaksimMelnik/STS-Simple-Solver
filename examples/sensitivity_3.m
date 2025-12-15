clearvars; close all; clc;
addpath('../data'); 

num_runs_sens = 60; 



%% выбор лучшего LHS
fprintf('=== Поиск лучшего LHS для устойчивой метамодели ===\n');

num_candidates = 2000;  % сколько разных LHS генерировать

d = 3;                 
p = 3;                
Nt = 20;                
T_steps = 501;               
n = num_runs_sens;

lP_0 = @(x) 1;
lP_1 = @(x) sqrt(3) * x;
lP_2 = @(x) sqrt(5)/2 * (3*x.^2 - 1); 
lP_3 = @(x) sqrt(7)/2 * (5*x.^3 - 3*x); 

best_score   = -Inf;
best_LHS  = [];
best_det  = -Inf;
best_cond = Inf;

% === 2. ПОИСК ЛУЧШЕГО LHS ===
for k = 1:num_candidates
    X = lhsdesign(n, d, 'criterion','maximin'); 
    X = 2 * X - 1;   % масштабируем в [-1,1]
  
    % СЛУ матрица Psi (NxM)
    Psi = zeros(n, Nt);
    
    for i = 1:n
        xi = X(i,:);
        
    % (20 полиномов)
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
    Psi(i, 11) = lP_3(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));    
    Psi(i, 12) = lP_0(xi(1)) * lP_3(xi(2)) * lP_0(xi(3));    
    Psi(i, 13) = lP_0(xi(1)) * lP_0(xi(2)) * lP_3(xi(3));    
    Psi(i, 14) = lP_2(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));    
    Psi(i, 15) = lP_2(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));    
    Psi(i, 16) = lP_1(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));    
    Psi(i, 17) = lP_0(xi(1)) * lP_2(xi(2)) * lP_1(xi(3));    
    Psi(i, 18) = lP_1(xi(1)) * lP_0(xi(2)) * lP_2(xi(3));    
    Psi(i, 19) = lP_0(xi(1)) * lP_1(xi(2)) * lP_2(xi(3));    
    Psi(i, 20) = lP_1(xi(1)) * lP_1(xi(2)) * lP_1(xi(3));    
     end
    
    G = Psi' * Psi;
    detG  = det(G);
    condG = cond(G);
    
    if condG < 10
        score = detG;      
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
results_matrix = zeros(T_steps, num_runs_sens + 1);  
 
for run_idx = 1:num_runs_sens
    fprintf('\n=== ЗАПУСК %d из %d ===\n', run_idx, num_runs_sens);
    
    current_point = LHS_points(run_idx, :);
    factors = 10.^(2 * current_point);
    
    load('../data/particles.mat', 'O2', 'O', 'Ar');
    
    O2.diss_Arrhenius_A('O')  = O2.diss_Arrhenius_A('O')  * factors(1);
    O2.diss_Arrhenius_A('O2') = O2.diss_Arrhenius_A('O2') * factors(2);
    O2.diss_Arrhenius_A('Ar') = O2.diss_Arrhenius_A('Ar') * factors(3);
    save 'particles1.mat' O2 O Ar;
    
    fprintf('Параметры запуска %d:\n', run_idx);
    fprintf('  O2+O:  %e\n', O2.diss_Arrhenius_A('O'));
    fprintf('  O2+O2: %e\n', O2.diss_Arrhenius_A('O2'));
    fprintf('  O2+Ar: %e\n', O2.diss_Arrhenius_A('Ar'));
    
    save('temp_sensitivity_vars.mat', 'run_idx', 'results_matrix', 'num_runs_sens', 'LHS_points', 'd', 'Nt', 'T_steps');
    
    fprintf('Запускаем O2_Ar_SW_Streicher...\n');
    O2_Ar_SW_Streicher;
    
    load('temp_sensitivity_vars.mat'); 
   
    lP_0 = @(x) 1; lP_1 = @(x) sqrt(3) * x;
    lP_2 = @(x) sqrt(5)/2 * (3*x.^2 - 1);
    lP_3 = @(x) sqrt(7)/2 * (5*x.^3 - 3*x); 
    
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


T_steps = size(results_matrix, 1);
time = results_matrix(:, 1);

%% Решение СЛУ  и Расчет Коэффициентов
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
    Psi(i, 11) = lP_3(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));    
    Psi(i, 12) = lP_0(xi(1)) * lP_3(xi(2)) * lP_0(xi(3));    
    Psi(i, 13) = lP_0(xi(1)) * lP_0(xi(2)) * lP_3(xi(3));    
    Psi(i, 14) = lP_2(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));    
    Psi(i, 15) = lP_2(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));    
    Psi(i, 16) = lP_1(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));    
    Psi(i, 17) = lP_0(xi(1)) * lP_2(xi(2)) * lP_1(xi(3));    
    Psi(i, 18) = lP_1(xi(1)) * lP_0(xi(2)) * lP_2(xi(3));    
    Psi(i, 19) = lP_0(xi(1)) * lP_1(xi(2)) * lP_2(xi(3));    
    Psi(i, 20) = lP_1(xi(1)) * lP_1(xi(2)) * lP_1(xi(3));    
end
alpha = zeros(Nt, T_steps);
for i = 1:T_steps
    R = results_matrix(i, 2:end)'; 
    alpha(:,i) = pinv(Psi) * R; 
end
%% Расчет индексов Соболя
IN_Sobol_Total = zeros(500, d + 1);    
IN_Sobol_Total(:,1) = time(2:501);

% Индексы полиномов, которые зависят от параметра 
Indices_T1 = [2, 5, 8, 9, 11, 14, 15, 16, 18, 20]; % X1 (A_{O2+O})
Indices_T2 = [3, 6, 8, 10, 12, 14, 16, 17, 19, 20]; % X2 (A_{O2+O2})
Indices_T3 = [4, 7, 9, 10, 13, 15, 17, 18, 19, 20]; % X3 (A_{O2+Ar})


for i = 2:T_steps 
    a = alpha(:, i);
    D = sum(a(2:end).^2); % Общая дисперсия (без первого коэффициента a1)
    
     D_T1 = sum(a(Indices_T1).^2);
     ST1 = D_T1 / D; 
  
     D_T2 = sum(a(Indices_T2).^2);
     ST2 = D_T2 / D; 
        
     D_T3 = sum(a(Indices_T3).^2);
     ST3 = D_T3 / D;
    
    
    IN_Sobol_Total(i-1, 2:4) = [ST1, ST2, ST3];
end

%% 
param_names = {'A_{O_2+O}', 'A_{O_2+O_2}', 'A_{O_2+Ar}'}; 
figure;
plot(IN_Sobol_Total(:,1), IN_Sobol_Total(:,2:4));
legend('A_O', 'A_O2', 'A_Ar');
xlabel('t');
ylabel('Total Sobol Indices');
%% ОЦЕНКА КАЧЕСТВА PCE (на обучающих точках)
T_steps = size(results_matrix, 1); 
all_rmse = zeros(1, T_steps);
all_nrmse = NaN(1, T_steps); 
sigma_threshold = 1e-10; 

fprintf('\n=== АНАЛИЗ КАЧЕСТВА PCE НА ОБУЧАЮЩЕМ НАБОРЕ ===\n');
for i = 1:T_steps
    R_true = results_matrix(i, 2:end)'; 
    R_PCE  = Psi * alpha(:, i); 
    
    rmse_t = sqrt(mean((R_true - R_PCE).^2));
    all_rmse(i) = rmse_t;
    
    sigma_Y_t = std(R_true); 
    
    if sigma_Y_t > sigma_threshold 
        nrmse_t = rmse_t / sigma_Y_t;
        all_nrmse(i) = nrmse_t; 
    end
end
% --- АНАЛИЗ И ВЫВОД РЕЗУЛЬТАТОВ ---
mean_nrmse = nanmean(all_nrmse);
fprintf('--------------------------------------------------------------------\n');
[max_nrmse, ~] = nanmax(all_nrmse);

fprintf('1.\n');
fprintf('   Макс. NRMSE: %.4f (%.2f%%)\n', max_nrmse, max_nrmse * 100);
fprintf('2.\n');
fprintf('   Среднее NRMSE: %.4f (%.2f%%) \n', mean_nrmse, mean_nrmse * 100);
fprintf('--------------------------------------------------------------------\n');

% 5. качественная оценка модели
if max_nrmse <= 0.20 && mean_nrmse <= 0.10
    fprintf(' КАЧЕСТВО МОДЕЛИ: Отличное \n');
elseif max_nrmse < 0.25 && mean_nrmse <= 0.15
    fprintf(' КАЧЕСТВО МОДЕЛИ: Приемлемое  \n');
else
    fprintf('КАЧЕСТВО МОДЕЛИ: Недостаточное \n');
     return;
end
   

%% ОЦЕНКА ОБОБЩАЮЩЕЙ СПОСОБНОСТИ НА ТЕСТОВОМ НАБОРЕ (Валидация)
% -------------------------------------------------------------------------
N_test = 200; % Количество тестовых запусков
% -------------------------------------------------------------------------
T_steps = size(results_matrix, 1);
time_val = zeros(T_steps, 1); 
R_true_tests = zeros(T_steps, N_test); 
R_PCE_tests  = zeros(T_steps, N_test); 

fprintf('\n\n=== 3. ПРОВЕРКА PCE НА %d НОВЫХ ТОЧКАХ ===\n', N_test);

% 1. Генерируем N_test новых, независимых точек X_test
X_test_norm = lhsdesign(N_test, d, 'criterion','maximin', 'iterations', 10);
X_test = 2 * X_test_norm - 1; 
save('lhs_test_plan.mat', 'X_test');

% 2. Цикл запуска Черного Ящика для всех тестовых точек
for run_idx = 1:N_test
    
    % СОХРАНЕНИЕ И ВОССТАНОВЛЕНИЕ ПЕРЕМЕННЫХ
    run_idx_val = run_idx; 
    save('pce_validation_state.mat', 'alpha', 'Nt', 'T_steps', 'd', 'lP_0', 'lP_1', 'lP_2', 'lP_3', 'run_idx_val', 'N_test');
    load('pce_validation_state.mat'); 
    load('lhs_test_plan.mat'); 

    fprintf('--> Запуск тестовой точки %d из %d...\n', run_idx_val, N_test); 
    current_point = X_test(run_idx_val, :); 
    factors = 10.^(2 * current_point);
    load('../data/particles.mat', 'O2', 'O', 'Ar');
    
    O2.diss_Arrhenius_A('O')  = O2.diss_Arrhenius_A('O')  * factors(1);
    O2.diss_Arrhenius_A('O2') = O2.diss_Arrhenius_A('O2') * factors(2);
    O2.diss_Arrhenius_A('Ar') = O2.diss_Arrhenius_A('Ar') * factors(3);
    save 'particles1.mat' O2 O Ar;
    
    O2_Ar_SW_Streicher;
    
    load('pce_validation_state.mat'); 
    load('lhs_test_plan.mat'); 
    
    % 3. Расчет PCE-предсказания
    time_val_single = time_ms_1; 
    R_true_single = n_O2_1 / Na; 
    
    xi = X_test(run_idx_val,:); 
    Psi_test_point = zeros(1, Nt); 
    
    % Расчет 20 базисных полиномов
    Psi_test_point(1, 1) = lP_0(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));      
    Psi_test_point(1, 2) = lP_1(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi_test_point(1, 3) = lP_0(xi(1)) * lP_1(xi(2)) * lP_0(xi(3)); 
    Psi_test_point(1, 4) = lP_0(xi(1)) * lP_0(xi(2)) * lP_1(xi(3)); 
    Psi_test_point(1, 5) = lP_2(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));  
    Psi_test_point(1, 6) = lP_0(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));  
    Psi_test_point(1, 7) = lP_0(xi(1)) * lP_0(xi(2)) * lP_2(xi(3)); 
    Psi_test_point(1, 8) = lP_1(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));  
    Psi_test_point(1, 9) = lP_1(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));  
    Psi_test_point(1, 10) = lP_0(xi(1)) * lP_1(xi(2)) * lP_1(xi(3));
    Psi_test_point(1, 11) = lP_3(xi(1)) * lP_0(xi(2)) * lP_0(xi(3));    
    Psi_test_point(1, 12) = lP_0(xi(1)) * lP_3(xi(2)) * lP_0(xi(3));    
    Psi_test_point(1, 13) = lP_0(xi(1)) * lP_0(xi(2)) * lP_3(xi(3));    
    Psi_test_point(1, 14) = lP_2(xi(1)) * lP_1(xi(2)) * lP_0(xi(3));    
    Psi_test_point(1, 15) = lP_2(xi(1)) * lP_0(xi(2)) * lP_1(xi(3));    
    Psi_test_point(1, 16) = lP_1(xi(1)) * lP_2(xi(2)) * lP_0(xi(3));    
    Psi_test_point(1, 17) = lP_0(xi(1)) * lP_2(xi(2)) * lP_1(xi(3));    
    Psi_test_point(1, 18) = lP_1(xi(1)) * lP_0(xi(2)) * lP_2(xi(3));    
    Psi_test_point(1, 19) = lP_0(xi(1)) * lP_1(xi(2)) * lP_2(xi(3));    
    Psi_test_point(1, 20) = lP_1(xi(1)) * lP_1(xi(2)) * lP_1(xi(3)); 
    
    R_PCE_single = (Psi_test_point * alpha)'; 
    
    temp_file_name = sprintf('temp_val_res_%d.mat', run_idx_val);
    save(temp_file_name, 'R_true_single', 'R_PCE_single', 'time_val_single');
    
end 

for run_idx = 1:N_test
    temp_file_name = sprintf('temp_val_res_%d.mat', run_idx);
    load(temp_file_name);
    R_true_tests(:, run_idx) = R_true_single;
    R_PCE_tests(:, run_idx) = R_PCE_single;
    time_val(:, 1) = time_val_single; 
    delete(temp_file_name); 
end
delete('pce_validation_state.mat'); 
delete('lhs_test_plan.mat'); 

% 5. Расчет ОБЩИХ метрик для ВСЕГО тестового набора
R_true_all = R_true_tests(:);
R_PCE_all = R_PCE_tests(:);
RMSE_test_total = sqrt(mean((R_true_all - R_PCE_all).^2));
sigma_Y_test_total = std(R_true_all);

if sigma_Y_test_total > 1e-10
    NRMSE_test_total = RMSE_test_total / sigma_Y_test_total;
else
    NRMSE_test_total = NaN; 
end

fprintf('\n--- ОБЩИЕ РЕЗУЛЬТАТЫ ВАЛИДАЦИИ (N=%d) ---\n', N_test);
fprintf('Общая NRMSE (на %d рядах): %.4f (%.2f%%)\n', N_test, NRMSE_test_total, NRMSE_test_total * 100);
fprintf('--------------------------------------------------------------------\n');

% 6. Оценка обобщения
if NRMSE_test_total < 0.10
    fprintf(' ОБОБЩЕНИЕ: Отличное.\n');
elseif NRMSE_test_total < 0.20
    fprintf(' ОБОБЩЕНИЕ: Приемлемое.\n');
else
    fprintf(' ОБОБЩЕНИЕ: Неудовлетворительное.\n');
end

%% Построение графика валидации
figure;

for run_idx = 1:N_test
    plot(time_val, R_true_tests(:, run_idx), '-', 'Color', colors(run_idx,:), 'LineWidth', 1.5);
    hold on;
    plot(time_val, R_PCE_tests(:, run_idx), '--', 'Color', colors(run_idx,:), 'LineWidth', 1);
end

xlabel('Время (ms)');
ylabel('Концентрация O2 (моль/м^3)');
title(sprintf('Валидация PCE: %d Тестовых Запусков (Общая NRMSE=%.2f%%)', N_test, NRMSE_test_total * 100));

h_legend = [plot(NaN, NaN, 'k-'), plot(NaN, NaN, 'k--')];
legend(h_legend, {'Истинный Черный Ящик', 'Предсказание PCE'}, 'Location', 'Best');
grid on;
hold off;