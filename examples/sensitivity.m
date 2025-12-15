clearvars; close all; clc;
addpath('../data');

% Количество запусков модели (вообще 10, но можно поэкспериментировать)
num_runs_sens = 30;


%% выбор лучшего LHS

num_candidates = 2000;   % сколько разных LHS генерировать
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
    
    %  1. Генерируем LHS
    X = lhsdesign(n, d, 'criterion','maximin', 'iterations', 10);  % в [0,1]
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

results_matrix = zeros(501, num_runs_sens + 1);  

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



%% Оценка качества PCE на обучающих точках (осреднение по времени)
T_steps = 501; 
time = results_matrix(:, 1);
all_rmse = zeros(1, T_steps);
all_nrmse = NaN(1, T_steps); 

% Задаем порог для стандартного отклонения (меньше него считаем артефактом)
sigma = 1e-10; 

fprintf('\n=== АНАЛИЗ КАЧЕСТВА PCE  ===\n');

for i = 1:T_steps
    % Истинные результаты для текущего временного шага
    R_true = results_matrix(i, 2:end)'; 
    
    % Предсказания PCE-модели
    R_PCE  = Psi * alpha(:, i); 
    
    % 1. Расчет абсолютной RMSE
    rmse_t = sqrt(mean((R_true - R_PCE).^2));
    all_rmse(i) = rmse_t;
    
    % 2. Расчет sigma_Y (стандартное отклонение)
    sigma_Y_t = std(R_true); 
    
    % 3. Расчет NRMSE с фильтрацией 
    if sigma_Y_t > sigma 
        nrmse_t = rmse_t / sigma_Y_t;
        all_nrmse(i) = nrmse_t; 
    end
end

% --- АНАЛИЗ И ВЫВОД РЕЗУЛЬТАТОВ ---

% Расчет среднего NRMSE, игнорируя NaN (то есть, t=0)
mean_nrmse = nanmean(all_nrmse);

fprintf(' Расчет качества PCE завершен для всех %d временных шагов.\n', T_steps);
% 1. max RMSE
[overall_max_rmse, overall_max_rmse_idx] = max(all_rmse);
time_at_overall_max_rmse = time(overall_max_rmse_idx);
nrmse_at_rmse_peak = all_nrmse(overall_max_rmse_idx);

fprintf('1. \n');
fprintf('   Макс. RMSE: %.6e\n', overall_max_rmse);
fprintf('   Достигнута при t = %.2f ms (Соответствующая NRMSE: %.4f)\n', time_at_overall_max_rmse, nrmse_at_rmse_peak);

% 2.max NRMSE
[max_nrmse, max_nrmse_idx] = nanmax(all_nrmse);
time_at_max_nrmse = time(max_nrmse_idx);
rmse_at_nrmse_peak = all_rmse(max_nrmse_idx); 

fprintf('2.\n');
fprintf('   Макс. NRMSE: %.4f (%.2f%%)\n', max_nrmse, max_nrmse * 100);
fprintf('   Достигнута при t = %.2f ms (Соответствующая RMSE: %.6e)\n', time_at_max_nrmse, rmse_at_nrmse_peak);


% 3. среднее значение (по времени) NRMSE
fprintf('3. \n');
fprintf('   Среднее NRMSE: %.4f (%.2f%%) (Исключая при t=0)\n', mean_nrmse, mean_nrmse * 100);

% 4. оценка при t = 0.0
fprintf('4. Оценка при t = 0.0 \n');
fprintf('   Абс. RMSE при t=0: %.6e\n', all_rmse(1));

% 5. качественная оценка модели
if max_nrmse <= 0.25 && mean_nrmse <= 0.20
    fprintf(' КАЧЕСТВО МОДЕЛИ: Отличное \n');
elseif max_nrmse < 0.30 && mean_nrmse <= 0.25 
    fprintf(' КАЧЕСТВО МОДЕЛИ: Приемлемое  \n');
else
    fprintf('КАЧЕСТВО МОДЕЛИ: Недостаточное \n');
     return;
    
end
%%  Оценка качества PCE на тестовом наборе (осреденение по всему пространству)

N_test = 100; % количество тестовых запусков
d = 3;
T_steps = 501; 
time_val = zeros(T_steps, 1); 
R_true_tests = zeros(T_steps, N_test); 
R_PCE_tests  = zeros(T_steps, N_test); 


fprintf('\n\n===  ПРОВЕРКА  НА PCE %d НОВЫХ ТОЧКАХ ===\n', N_test);

% 1. Генерируем N_test новых, независимых точек X_test
X_test_norm = lhsdesign(N_test, d, 'criterion','maximin', 'iterations', 10);
X_test = 2 * X_test_norm - 1; 

% 2. Цикл запуска Черного Ящика для всех тестовых точек
for run_idx = 1:N_test
    
    
    run_idx_val = run_idx; 
    
    save('pce_validation_state.mat', 'alpha', 'Nt', 'X_test', 'T_steps', 'd', 'lP_0', 'lP_1', 'lP_2', 'run_idx_val', 'N_test','LHS_points');

    load('pce_validation_state.mat'); 

    fprintf('--> Запуск тестовой точки %d из %d...\n', run_idx_val, N_test); 
    
    current_point = X_test(run_idx_val, :); 
    factors = 10.^(2 * current_point);

    % Загружаем и меняем исходные данные для запуска
    load('../data/particles.mat', 'O2', 'O', 'Ar');
    
    O2.diss_Arrhenius_A('O')  = O2.diss_Arrhenius_A('O')  * factors(1);
    O2.diss_Arrhenius_A('O2') = O2.diss_Arrhenius_A('O2') * factors(2);
    O2.diss_Arrhenius_A('Ar') = O2.diss_Arrhenius_A('Ar') * factors(3);
    save 'particles1.mat' O2 O Ar;
   
    O2_Ar_SW_Streicher;
    
    load('pce_validation_state.mat'); 
    
    
    % Сохраняем результаты 
    time_val_single = time_ms_1; 
    R_true_single = n_O2_1 / Na; 
    
    % Расчет PCE-предсказания
    xi = X_test(run_idx_val,:); 
    Psi_test_point = zeros(1, Nt); 
    
    % Расчет базисных полиномов для текущей тестовой точки
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
    
    R_PCE_single = (Psi_test_point * alpha)'; 
    
    % Сохраняем результаты этого шага в отдельный файл
    temp_file_name = sprintf('temp_val_res_%d.mat', run_idx_val);
    save(temp_file_name, 'R_true_single', 'R_PCE_single', 'time_val_single');
    
end

% 4.заполняем результатами
for run_idx = 1:N_test
    temp_file_name = sprintf('temp_val_res_%d.mat', run_idx);
    load(temp_file_name);
    R_true_tests(:, run_idx) = R_true_single;
    R_PCE_tests(:, run_idx) = R_PCE_single;
    time_val(:, 1) = time_val_single; % Временный вектор
    delete(temp_file_name); % Удаляем временный файл
end
delete('pce_validation_state.mat'); % Удаляем файл 

% 5. Расчет ОБЩИХ метрик 
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
fprintf('Общая RMSE (на %d рядах): %.6e\n', N_test, RMSE_test_total);
fprintf('Общая NRMSE (на %d рядах): %.4f (%.2f%%)\n', N_test, NRMSE_test_total, NRMSE_test_total * 100);

% 6. Оценка обобщения
if NRMSE_test_total < 0.10
    fprintf('ОБОБЩЕНИЕ: Отличное.\n');
elseif NRMSE_test_total < 0.20
    fprintf('ОБОБЩЕНИЕ: Приемлемое. \n');
else
    fprintf('ОБОБЩЕНИЕ: Неудовлетворительное. \n');
end

% 7. Построение графика
figure;
colors = jet(N_test); % Цветовая карта для разных запусков

for run_idx = 1:N_test
    % Истинный ряд (сплошная линия)
    plot(time_val, R_true_tests(:, run_idx), '-', 'Color', colors(run_idx,:), 'LineWidth', 1.5);
    hold on;
    % Предсказание PCE (пунктирная линия того же цвета)
    plot(time_val, R_PCE_tests(:, run_idx), '--', 'Color', colors(run_idx,:), 'LineWidth', 1);
end

xlabel('Время (ms)');
ylabel('Концентрация O2 (моль/м^3)');
title(sprintf('Валидация PCE: %d Тестовых Запусков (Общая NRMSE=%.2f%%)', N_test, NRMSE_test_total * 100));

h_legend = [plot(NaN, NaN, 'k-'), plot(NaN, NaN, 'k--')];
legend(h_legend, {'Истинный Черный Ящик', 'Предсказание PCE'}, 'Location', 'Best');

grid on;
hold off;