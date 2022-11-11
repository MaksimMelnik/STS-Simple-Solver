% код по загрузке промежуточных данных, сохранённых по способу ОВ
% 17.12.2020

switch ind
    case 1
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_1Torr_O2X_1\ChemDensitiesTime.txt');
    case 2
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_1Torr_O2X_0p7\ChemDensitiesTime.txt');
    case 3
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_10Torr_O2X_1\ChemDensitiesTime.txt');
    case 4
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_10Torr_O2X_0p7\ChemDensitiesTime.txt');
end
load print_data_test.txt
tm = print_data_test(:,1); % MP - массив загруженных данных, tm - вектор с шагами интегрирования

tm_i = []; % массив с значащими индексами
for i = 2:length(tm)-1
    ti = tm(i); tip = tm(i+1);
    if ti == tip
        continue;
    else
        tm_i = [tm_i; i];
        while tm(tm_i(end))>tip
            tm_i=tm_i(1:end-1);
        end
    end
end

win_1=2; win_2=3;
figure('Position', [70 0 1200 610])
subplot(win_1, win_2, 1)
semilogx(data3(:,1), data3(:,2), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O2), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])%6e-2])%tm(end)])
title('nO2(x)')
subplot(win_1, win_2, 2)
loglog(data3(:,1), data3(:,3), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O2a), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO2(a)')
subplot(win_1, win_2, 3)
loglog(data3(:,1), data3(:,4), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O2b), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO2(b)')
subplot(win_1, win_2, 4)
loglog(data3(:,1), data3(:,6), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO(3P)')
subplot(win_1, win_2, 5)
loglog(data3(:,1), data3(:,7), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O1D), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO(1D)')
subplot(win_1, win_2, 6)
loglog(data3(:,1), data3(:,10), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O3), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO3')

figure('Position', [70 0 1200 610])
win_1=2; win_2=2;
subplot(win_1, win_2, 1)
loglog(data3(:,1), data3(:,5), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O2p), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO2+')
subplot(win_1, win_2, 2)
loglog(data3(:,1), data3(:,8), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_Op), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO+')
subplot(win_1, win_2, 3)
loglog(data3(:,1), data3(:,9), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_Om), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO-')
subplot(win_1, win_2, 4)
loglog(data3(:,1), data3(:,11), 'LineWidth', 1.5)
hold on 
loglog(tm(tm_i), print_data_test(tm_i,1+num_O3exc), 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 tm(end)])
title('nO3(exc)')